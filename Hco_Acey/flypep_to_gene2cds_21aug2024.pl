#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my @outputs = ();

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) (.*) \z/xms ) { 
        my $protein = $1;
        my $annot   = $2;
        if ( $annot =~ / name=(\S+?)\-P[A-Z] /xms ) {
            my $gene = $1;
            my $output = "$gene\t$protein";
            push @outputs, $output;
            $output = q{};
        }
        if ( $annot =~ / FlyBase_Annotation_IDs:(\S+?)\-P[A-Z] /xms ) {
            my $gene = $1;
            my $output = "$gene\t$protein";
            push @outputs, $output;
            $output = q{};
        }
    }
}

@outputs = sort(@outputs);
@outputs = uniq(@outputs);

foreach my $output (@outputs) {
    if ( $output =~ /\S/xms ) {
        print "$output\n";
    }
}
