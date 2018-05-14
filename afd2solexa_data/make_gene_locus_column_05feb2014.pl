#!/usr/bin/env perl

use strict;
use warnings;

my $header = "Gene\tMap_position";
my @outputs = ();

while (my $input = <>) { 
    chomp $input;
    $input =~ s/["]//g;
    if ( $input =~ /\A (WBGene\d+) \t (\S+) \t (\S+) \t (\S+) /xms ) { 
        my $gene     = $1;
        my $chr      = $2;
        my $map      = $3;
        my $err      = $4;
        my $out_text = $gene . "\t" . $chr . q{: } . $map . q{ [+/-] } . $err;
        push @outputs, $out_text;
    }
}

foreach my $output (@outputs) { 
    print "$header\n" if $header;
    $header = q{};
    print "$output\n";
}

