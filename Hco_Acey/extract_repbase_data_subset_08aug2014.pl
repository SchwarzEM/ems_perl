#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $list = $ARGV[0];
my $data = $ARGV[1];

my %wanted = ();

my $header = "Element name\tSpecies\tElement type\tSource";

open my $LIST, '<', $list;
while (my $input = <$LIST>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) { 
        $wanted{$input} = 1;
    }
    else { 
        die "From list file $list, can't parse input: $input\n";
    }
}
close $LIST;

open my $DATA, '<', $data;
while (my $input = <$DATA>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \b/xms ) {
        my $seq = $1;
        if ( exists $wanted{$seq} ) {
            print "$header\n" if $header;
            $header = q{};
            print "$input\n";
        }
    }
    else {
        die "From data file $data, can't parse input: $input\n";
    }
}
close $DATA;

