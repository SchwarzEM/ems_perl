#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A InterPro \s+ (IPR\d+) \s+ (\S.+\S) \s* \z/xms ) { 
        my $desc = $2;
        my $ipr  = $1;
        my $full_desc = "$desc [$ipr]";
        $data_ref->{'interpro'}->{$full_desc} = 1;
    }
}

my @interpros = sort keys %{ $data_ref->{'interpro'} };

foreach my $interpro1 (@interpros) { 
    print "$interpro1\n";
}

