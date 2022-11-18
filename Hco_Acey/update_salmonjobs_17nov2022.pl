#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $infile_list = q{};
$infile_list    = $ARGV[0] if $ARGV[0];

if (! -e $infile_list ) {
    die "Format: update_salmonjobs_17nov2022.pl [list of input scripts] -- ",
        "script prints out list of renamed/updated scripts in working directory\n";
}

my @infiles = ();

open my $LIST, '<', $infile_list;
while (my $infile = <$LIST>) {
    chomp $infile;
    if (! -e $infile) {
        die "Nonexistent input file: $infile\n"
    }
    my $outfile = $infile;
    $outfile    = basename($outfile);
    $outfile    =~ s/2022\.02\.12/2022.11.18/;
    $outfile    = safename($outfile);

    open my $INFILE, '<', $infile;
    open my $OUTFILE, '>', $outfile;
    while (my $input = <$INFILE>) {
        chomp $input;
        $input =~ s/2022\.02\.12\./2022.11.18./g;
        $input =~ s/salmon_1\.6\.0/salmon_1.9.0/g;
        $input =~ s/Acey\.v2_WBPS16_gentrome_index/Acey.v2.1_WBPS16_gentrome_index/g;
        $input =~ s/ \-\-validateMappings//g;
        $input =~ s#Acey/2022\.02\.03/annots/Acey_v2\.2022\.02\.02\.01\.cds2gene\.tsv\.txt#Acey/2022.02.13/annots/Acey_v2.1.2022.11.14.01.cds2gene.tsv.txt#g;
        print $OUTFILE "$input\n";        
    }
}
close $LIST;

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

