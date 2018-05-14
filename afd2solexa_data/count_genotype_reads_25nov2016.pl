#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (.+) [ ] rep \. [ ] \d \t (\S+) \t/xms ) {
        my $genotype = $1;
        my $counts   = $2;
        $counts =~ s/[,]//g;
        $data_ref->{'genotype'}->{$genotype}->{'counts'} += $counts;
    }
}

my @genotypes = sort keys %{ $data_ref->{'genotype'} };
foreach my $genotype (@genotypes) {
    my $counts = $data_ref->{'genotype'}->{$genotype}->{'counts'};
    $counts    = commify($counts);
    print "$genotype\t$counts\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

