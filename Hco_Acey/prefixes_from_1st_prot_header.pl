#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

my $proteome = q{};
my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ==> \s+ (\S+) \s+ /xms ) { 
        $proteome = $1;
        $proteome = basename($proteome);
    }
    elsif ( $proteome and ( $input =~ /\A > ([A-Za-z0-9]+) _ /xms ) ) { 
        my $prefix = $1;
        if (     ( exists $data_ref->{'prefix'}->{$prefix}->{'proteome'} ) 
             and ( $data_ref->{'prefix'}->{$prefix}->{'proteome'} ne $proteome ) ) {
            die "Inconsistent mapping of prefix $prefix to two proteomes: $proteome vs. $data_ref->{'prefix'}->{$prefix}->{'proteome'}\n";
        }
        $data_ref->{'prefix'}->{$prefix}->{'proteome'} = $proteome;
    }
}

my @prefixes = sort keys %{ $data_ref->{'prefix'} };
foreach my $prefix (@prefixes) { 
    my $proteome1 = $data_ref->{'prefix'}->{$prefix}->{'proteome'};
    print "    '$prefix' => '$proteome1',\n";
}


