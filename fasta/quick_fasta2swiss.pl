#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Bio::SeqIO;

foreach (my $infile = <>) { 
    chomp $infile;
    my $outfile = $infile ;
    $outfile =~ s/\.fa\z//;
    $outfile = "$outfile.swiss";
    $outfile = safename($outfile);
    my $infile_object  = Bio::SeqIO->new(-file => "$infile",   -format => 'fasta');
    my $outfile_object = Bio::SeqIO->new(-file => ">$outfile", -format => 'swiss');
    while ( my $seq_data = $infile_object->next_seq() ) {
        $outfile_object->write_seq($seq_data);
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

