#!/usr/bin/perl

# gffdna2prot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/25/2008.
# Purpose: from Bio:DB:GFF with *.fa/*.gff of pred. CDSes, get proteome.

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

# Use a simple sort; not WashU-aware.

my @contigs = sort 
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

