#!/usr/bin/perl

# get_Etch_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/24/2007.
# Purpose: get Etchberger subset of genes, then use it to filter Ali's data.

use strict;
use warnings;

unless ($#ARGV == 1 ) { 
    die "Format: ./get_Etch_subset.pl etch.csv ali_data.txt\n";
}

my $etch     = $ARGV[0];
my $ali_data = $ARGV[1];
my %etch_cds = ();

open (my $ETCH, "<", $etch) or 
    die "Can't open $etch: $!";
while (my $input = <$ETCH>) { 
    chomp $input;
    if ( $input =~ / \A [^\"]+ \" ( [A-Z0-9]+\.[A-Z0-9]+ ) \" /xms ) {
        $etch_cds{$1} = 1;
    }
}
close $ETCH;

open (my $ALI_DATA, "<", $ali_data) or 
    die "Can't open $ali_data: $!";
while (my $input = <$ALI_DATA>) { 
    chomp $input;
    if ( ( $input =~ / \A WBGene\d+  [\|\s]+ ( [A-Z0-9]+\.[A-Z0-9]+ ) /xms ) 
          and ( $etch_cds{$1} ) ) { 
        print "$input\n";
    }
}
close $ALI_DATA;

