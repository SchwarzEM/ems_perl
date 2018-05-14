#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @infiles = @ARGV;

foreach my $infile (@infiles) {
    chomp $infile;
    if ( $infile =~ /\A (\S+) \.tsv\.txt \z/xms ) {
        my $stem    = $1;
        my $outfile = "$stem.subset.tsv.txt";
        $outfile    = safename($outfile);
        open my $INFILE, '<', $infile;
        open my $OUTFILE, '>', $outfile;
        while (my $input = <$INFILE>) {
            chomp $input;
            if ( $input =~ /\A ([^\t]*) \t ([^\t]*) \t [^\t]* \t [^\t]* \t ([^\t]*) \z/xms ) {
                my $gene   = $1;
                my $log2fc = $2;
                my $fdr    = $3;
                print $OUTFILE "$gene\t$log2fc\t$fdr\n";
            }
            else {
                die "From input file $infile, cannot parse: $input\n";
            }
        }
        close $INFILE;
        close $OUTFILE;
    }
    else {
        die "Cannot parse input file name: $infile\n";
    }
}

# Purpose: always print/export data to a filename that's new, avoiding overwriting existing files.
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

