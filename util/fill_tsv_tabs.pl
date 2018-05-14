#!/usr/bin/env perl

# fill_tsv_tabs.pl -- Erich Schwarz <emsch@@its.caltech.edu>, 11/23/2010.
# Purpose: given a ragged, uneven *.tsv produced by Microsoft's Excel, scan for the true tab-count and then print a version with full tabs per line.

use strict;
use warnings;

my @input_files = @ARGV;

foreach my $infile (@input_files) { 
    my $tabcount = 0;
    open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
    while ( my $input = <$INFILE>) { 
        my $line_tabcount = ( $input =~ tr/\t/\t/ );
        if ( $line_tabcount > $tabcount ) { 
            $tabcount = $line_tabcount;
        }
    }
    close $INFILE or die "Can't close filehandle to $infile: $!";
    open $INFILE, '<', $infile or die "Can't open input file $infile: $!";
    while ( my $input = <$INFILE>) {
        chomp $input;
        my $line_tabcount = ( $input =~ tr/\t/\t/ );
        if ( $line_tabcount < $tabcount ) {
            my $extra_tabcount = ( $tabcount - $line_tabcount );
            my $added_tabs = ("\t" x $extra_tabcount);
            $input .= $added_tabs;
        }
        print "$input\n";
    }
    close $INFILE or die "Can't close filehandle to $infile: $!";
}

