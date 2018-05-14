#!/usr/bin/env perl

# gffdna2prot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/8/2008.
# Purpose: from Bio:DB:GFF with *.fa/*.gff of pred. CDSes, get proteome.
# Note: supports two different sortings of contig/chrom. names.

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::GFF;
use Getopt::Long;

my ($CODING_UNIT, $FEATURE_TYPE, $AGGREGATOR, $QUERY_BASE, $SORT_TYPE);

GetOptions ( 'coding_unit:s' => \$CODING_UNIT,
             'type:s'        => \$FEATURE_TYPE,
             'query_base:s'  => \$QUERY_BASE,
             'aggregator:s'  => \$AGGREGATOR, 
             'sort_type:s'   => \$SORT_TYPE, );

if (! $QUERY_BASE) {
    die "Minimum args.: --query_base/-q [queried database]\n";
}

# Local defaults for most options:

$CODING_UNIT  ||= 'CDS';
$FEATURE_TYPE ||= 'contig:assembly';
$AGGREGATOR   ||= 'processed_transcript';
$SORT_TYPE    ||= 'WashU-aware';

my $dna = q{};
my $db = Bio::DB::GFF->new( -dsn => $QUERY_BASE); 

my $prot_file = join('.', &get_local_date()) . ".prot.fa";
if (-e $prot_file) {
    die "Will not overwrite $prot_file\n";
}
my $seq_out 
    = Bio::SeqIO->new( -file => ">>$prot_file", -format => 'Fasta', );

my @contigs = ();

if ( $SORT_TYPE eq 'WashU-aware' ) { 
    @contigs = sort { &numWUCont($a) <=> &numWUCont($b) } 
               map { $_->display_id }
               $db->features( -types => $FEATURE_TYPE );
}
if ( $SORT_TYPE eq 'plain' ) {
    @contigs = sort 
               map { $_->display_id }
               $db->features( -types => $FEATURE_TYPE );
}
if ( ( $SORT_TYPE ne 'WashU-aware' )
     and ( $SORT_TYPE ne 'plain' ) ) {
    die "Sort type must be 'WashU-aware' or 'plain', not '$SORT_TYPE'!\n";
}

foreach my $contig (@contigs) { 
    my $segment1 = $db->segment($contig);
    my @p_txs = $segment1->features($AGGREGATOR);
    foreach my $p_t (sort @p_txs) {
        $dna = q{};
        my @CDSes = $p_t->$CODING_UNIT;

        my $cds_name = $CDSes[0]->display_id();

        # $cds->seq == Bio::PrimarySeq, *not* clean nt seq.!
        foreach my $cds (@CDSes) { 
            $dna = $dna . $cds->seq->seq;  
        }
        my $full_cds = Bio::Seq->new( -display_id => $cds_name, 
                                      -seq => $dna, );
        my $prot = $full_cds->translate;
        $seq_out->write_seq($prot);
    }
}

sub get_local_date {
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;
}

sub numWUCont { 
     my $contig = $_[0];
     if ( ( $contig !~ / Contig0*[1-9]\d* \z/xms ) 
           and   ( $contig !~ /\A Contig0 \z/xms ) ) {
         die "Malformatted input.\n";
     }
     if ($contig =~ /\A Contig0 \z/xms) {
         $contig = 0;
     }
     # Strip zeros, avoid octal-izing
     if ($contig =~ /\A Contig0*([1-9]\d*) \z/xms) { 
         $contig = $1;
     }
     return $contig;
}

