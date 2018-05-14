#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t .+ \t (\S+) \t \S+ \z/xms ) {
        my $query = $1;
        my $e_val = $2;
        if ( exists $data_ref->{'query'}->{$query}->{'e_val'} ) {
            my $prior_e = q{};
            $prior_e = $data_ref->{'query'}->{$query}->{'e_val'};
            if ( $prior_e > $e_val ) {
                $data_ref->{'query'}->{$query}->{'e_val'} = $e_val;
            }
        }
        else {
            $data_ref->{'query'}->{$query}->{'e_val'} = $e_val;
        }
    }
    elsif ( $input !~ /\A [#] /xms ) {
        die "Cannot parse input line: $input\n";
    }
}

my @queries = sort keys %{ $data_ref->{'query'} };
foreach my $query (@queries) {
    my $e_val = $data_ref->{'query'}->{$query}->{'e_val'};
    print "$query\t$e_val\n";
}


