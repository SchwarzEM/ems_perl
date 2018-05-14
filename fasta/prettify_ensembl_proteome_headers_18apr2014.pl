#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix = $ARGV[0];
my $infile = $ARGV[1];

if ( (! $prefix) or (! $infile) ) {
    die "Format: prettify_ensembl_proteome_headers_18apr2014.pl [prefix, usefully for species] [input proteome file, or stream '-'] > [outfile proteome]\n";
}

my $INFILE;

if ($infile eq '-') {
    # Special case: get the stdin handle
    $INFILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open my $INFILE, '<', $infile;
}

while (my $input = <$INFILE>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+ .* gene:(\S+)) /xms ) { 
            my $header = $1;
            my $gene   = $2;
            $gene = $prefix . $gene;
            print ">$gene  $header\n";
        }
        else {
            die "Can't parse header: $input\n";
        }
    }
    else { 
        print "$input\n";
    }
}

close $INFILE;

