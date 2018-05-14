#!/usr/bin/perl

# gffdna2prot_WUaware.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/17/2008.
# Purpose: from Bio:DB:GFF with *.fa/*.gff of pred. CDSes, get proteome, w/ WashU-aware sorting.

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::GFF;

my $query_base = $ARGV[0];
my $dna = q{};
my $db = Bio::DB::GFF->new( -dsn => $query_base); 

my $prot_file = join('.', &get_local_date()) . ".prot.fa";
if (-e $prot_file) {
    die "Will not overwrite $prot_file\n";
}
my $seq_out 
    = Bio::SeqIO->new( -file => ">>$prot_file", -format => 'Fasta', );

my @contigs = sort { &numWUCont($a) <=> &numWUCont($b) } 
              map { $_->display_id }
              $db->features( -types => 'contig:assembly' );

foreach my $contig (@contigs) { 
    my $segment1 = $db->segment($contig);
    my @p_txs = $segment1->features('processed_transcript');
    foreach my $p_t (sort @p_txs) {
        $dna = q{};
        my @CDSes = $p_t->CDS;

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

