#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @input_uniprots = @ARGV;

my $miss_list = pop @input_uniprots;

my %seen = ();

open my $MISS, '<', $miss_list;
while (my $input = <$MISS>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) { 
        $seen{$input} = 1;
    }
    else { 
        die "From miss-list $miss_list, can't parse: $input\n";
    }
}
close $MISS;

foreach my $uniprot (@input_uniprots) { 
    open my $UPROT, '<', $uniprot;
    while (my $input = <$UPROT>) {
        chomp $input;
        if ( $input =~ /\A [A-Z] [a-z]+ \t (\S+) \z/xms ) {
            my $protein = $1;
            if ( exists $seen{$protein} ) { 
                print "$input\n";
            }
        }
        else {
            die "From UniProt-based taxonomy-tag table $uniprot, can't parse: $input\n";
        }
    }
    close $UPROT;
}

