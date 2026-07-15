#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $es_list = q{};
my $mapping = q{};

$es_list = $ARGV[0] if $ARGV[0];
$mapping = $ARGV[1] if $ARGV[1];

my %esp = ();

if ( (! $es_list ) or (! $mapping ) ) {
    die "Format: filt_WashU.ES_15jul2026a.pl [ESP gene list] [mapping] > [subset of mapping with ESP genes]\n";
}

open my $ES, '<', $es_list;
while ( my $gene = <$ES> ) {
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        $esp{$gene} = 1;
    }
    else {
        die "From ESP gene list $es_list, cannot parse: $gene\n";
    }
}
close $ES;

open my $MAP, '<', $mapping;
while ( my $input = <$MAP> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \t (\S+) \z/xms ) {
        my $gene = $1;
        if ( exists $esp{$gene} ) {
            print "$input\n";
        }
    }
}
close $MAP;

