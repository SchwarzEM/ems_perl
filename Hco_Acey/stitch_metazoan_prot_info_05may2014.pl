#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $selection  = $ARGV[0];
my $data_table = $ARGV[1];

my %choices = (
    'Name' => 1,
);

open my $SELECT, '<', $selection;
while (my $input = <$SELECT>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) { 
        $choices{$input} = 1;
    }
    else {
        die "Can't parse input line from selection table $selection: $input\n";
    }
}
close $SELECT;

open my $DATA, '<', $data_table;
while (my $input = <$DATA>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s /xms ) {
        my $name = $1;
        if (exists $choices{$name} ) { 
            print "$input\n";
        }
    }
    else {
        die "Can't parse input line from selection table $selection: $input\n";
    }
}
close $DATA;

