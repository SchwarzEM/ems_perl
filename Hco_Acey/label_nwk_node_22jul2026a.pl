#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

if (! @ARGV ) {
    die "Usage: $0 [Newick file to be internally tagged] [taxon 1 to be internally tagged] [(optionally, 2-N more taxa)]\n";
}

my $newick = shift @ARGV;
my @taxa   = @ARGV;
@ARGV      = ($newick);

$^I = '.prev';

while ( my $input = <> ) {
    chomp $input;
    foreach my $taxon (@taxa) {
        $input =~ s/$taxon/$taxon\{$taxon\}/g;
    }
    print "$input\n";
}
