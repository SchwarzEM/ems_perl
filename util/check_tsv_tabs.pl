#!/usr/bin/env perl

# check_tsv_tabs.pl -- Erich Schwarz <emsch@@its.caltech.edu>, 11/23/2010.
# Purpose: check one or more *.tsv files for absolutely consistent tabcounts in each line of each file.

use strict;
use warnings;

my $tabcount = 0;
my @input_files = @ARGV;

my @good_files = ();
my @bad_files  = ();

LOOP: foreach my $infile (@input_files) { 
    open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
    $tabcount = 0;
    while ( my $input = <$INFILE>) { 
        my $line_tabcount = ( $input =~ tr/\t/\t/ );
        if ( $line_tabcount != $tabcount ) { 
            if ( $tabcount > 0 ) { 
                my $bad_file = "    $infile has inconsistent tabcount (first $tabcount, then $line_tabcount)\n";
                push @bad_files, $bad_file;
                next LOOP;
            }
            $tabcount = $line_tabcount;
        }
    }
    close $INFILE or die "Can't close filehandle to $infile: $!";
    my $good_file = "    $infile has consistent tabcount of $tabcount per line\n";
    push @good_files, $good_file;
}

my $good_file_count = @good_files;
my $bad_file_count  = @bad_files;

print "$good_file_count good files:\n";
print @good_files;
print "\n";
print "$bad_file_count bad files:\n";
print @bad_files;

