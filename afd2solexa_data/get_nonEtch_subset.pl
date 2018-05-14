#!/usr/bin/perl

# get_Etch_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/24/2007.
# Purpose: get Etchberger subset of genes, then use it to filter Ali's data.

use strict;
use warnings;

unless ($#ARGV == 2 ) { 
    die "Format: ./get_Etch_subset.pl  etch.csv  etch2.csv  ali_data.txt\n";
}

my @etches    = ();
$etches[0] = $ARGV[0];
$etches[1] = $ARGV[1];
my $ali_data  = $ARGV[2];
my %etch_cds  = ();

foreach my $etch (@etches) { 
    open (my $ETCH, "<", $etch) or 
        die "Can't open $etch: $!";
    while (my $input = <$ETCH>) { 
        chomp $input;
        if ( $input =~ / \A [^\"]+ \" ( [A-Z0-9]+\.[A-Z0-9]+ ) \" /xms ) {
            $etch_cds{$1} = 1;
        }
    }
    close $ETCH;
}

open (my $ALI_DATA, "<", $ali_data) or 
    die "Can't open $ali_data: $!";
while (my $input = <$ALI_DATA>) { 
    chomp $input;
    if ( ( $input =~ / \A WBGene\d+ \| (?: \S+ \| )? ( [A-Z0-9]+\.[A-Z0-9]+ ) /xms ) 
          and (! $etch_cds{$1} ) ) { 
        print "$input\n";
    }
}
close $ALI_DATA;

