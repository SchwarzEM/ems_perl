#!/usr/bin/env perl

# simple_fastq_splitter.pl -- Erich Schwarz <ems394@cornell.edu>, 8/2/2016.
# Purpose: very simple-minded splitter of interleaved paired-end FastQ file into two different paired-end FastQ files; does not deal with non-stereotypically formatted FastQ but should not jam on third-line variants.

use strict;
use warnings;
use autodie;

my $infile   = $ARGV[0];
my $outfile1 = $ARGV[1];
my $outfile2 = $ARGV[2];

$outfile1 = safename($outfile1);
$outfile2 = safename($outfile2);

open my $INFILE,  '<', $infile;
open my $OUTFILE, '>', $outfile1;

my $i = 0;

while (my $input = <$INFILE>) {
    chomp $input;
    $i++;
    if ( ($i % 8) == 1 ) {
        close $OUTFILE;
        # *append* more lines; don't *overwrite* them, with '>' in this next line:
        open  $OUTFILE, '>>', $outfile1;
        print $OUTFILE "$input\n";
    }
    elsif ( ($i % 8) == 5 ) {
        close $OUTFILE; 
        # *append*, don't overwrite:
        open  $OUTFILE, '>>', $outfile2;
        print $OUTFILE "$input\n";
    }
    else {
        print $OUTFILE "$input\n";
    }
}

close $INFILE;
close $OUTFILE;

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

