#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix = $ARGV[0];
my $infile = $ARGV[1];

my $print = 0;

if ( (! $prefix) or ( $prefix !~ /\A \S+ \z/xms ) or (! $infile) ) {
    die "Format: extract_fasta_prefix.pl [seqname prefix to select; no spaces] [FASTA seq. input file] > [STDOUT: selected FASTA subset]\n";
}

open my $INFILE, '<',  $infile;
while (my $input = <$INFILE>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) {
        # No use of 'xms' here: 
        if ( $input =~ /\A>$prefix/ ) { 
            print "$input\n";
            $print = 1;
        }
        else {
            $print = 0;
        }
    }
    elsif ($print) { 
        print "$input\n";
    }
}
close $INFILE;

