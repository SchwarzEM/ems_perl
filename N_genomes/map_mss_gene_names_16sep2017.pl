#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $first_names = q{};
my $second_names = q{};

$first_names = $ARGV[0] if $ARGV[0];
$second_names = $ARGV[1] if $ARGV[1];

my %cds2gene = ();
my %old2new = ();

open my $FIRST, '<', $first_names;
while (my $input = <$FIRST>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) { 
        my $old_gene = $1;
        my $cds      = $2;
        $cds2gene{$cds} = $old_gene;
    }
    else {
        die "In first file $first_names, cannot parse: $input\n";
    }
}
close $FIRST;

open my $SECOND, '<', $second_names;
while (my $input = <$SECOND>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $cds      = $1;
        my $new_gene = $2;
        if ( exists $cds2gene{$cds} ) {
            my $old_gene = $cds2gene{$cds};
            $old2new{$old_gene} = $new_gene;
        }
        else {
            $old2new{$cds} = $new_gene;
        }
    }
    else {
        die "In second file $second_names, cannot parse: $input\n";
    }
}

my @old_names = sort keys %old2new;

foreach my $old_name (@old_names) {
    my $new_name = $old2new{$old_name};
    print "$old_name\t$new_name\n";
}


