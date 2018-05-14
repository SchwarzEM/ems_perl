#!/usr/bin/env perl

# split_ncbi_tbl.pl -- Erich Schwarz <ems@emstech.org>, 2/23/2014.
# Given an NCBI feature table meant for tbl2sqn, split it into individual tables, each one having the name of its intended sequence and the suffix ".tbl" (subject to safenaming).

use strict;
use warnings;
use autodie;
use Scalar::Util qw(openhandle);

my $sequence = q{};
my $outfile  = q{};

my $TABLE;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > Features \s+ (\S+) \b .* \z/xms ) { 
        # Detect any new sequence name.
        $sequence = $1;

        # Stop printing to whatever previous file was being printed to.
        close $TABLE if openhandle $TABLE;

        # Define a new file to print to, and open printing to it.
        $outfile = $sequence . ".tbl";
        $outfile = safename($outfile);
        open $TABLE, '>', $outfile;
        print $TABLE "$input\n";
    }
    else {
        print $TABLE "$input\n";
    }
}

# For the sake of completeness:
close $TABLE if openhandle $TABLE;

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


