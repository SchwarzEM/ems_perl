#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my @orig_annots = ();

while (my $input = <>) {
    chomp $input;
    push @orig_annots, $input;
}

my @sorted_annots = sort { annot_val($a) <=> annot_val($b) } @orig_annots;

foreach my $sorted_annot (@sorted_annots) {
    print "$sorted_annot\n";
}

sub annot_val {
    my $_input = $_[0];
    chomp $_input;
    if ( $_input =~ /\A [^\t]* \t (\S+) \t /xms ) {
        my $_annot_val = $1;
        if (! looks_like_number($_annot_val) ) {
            die "Putative annotation value \"$_annot_val\" is apparently non-numerical in: $_input\n";
        }
        return $_annot_val;
    }
    else {
        die "Cannot parse annotation value in: $_input\n";
    }
    return;
}

