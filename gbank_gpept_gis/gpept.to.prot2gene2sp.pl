#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gpept   = q{};
my $species = q{};

my $protein = q{};
my $gene    = q{};

$gpept   = $ARGV[0] if $ARGV[0];
$species = $ARGV[1] if $ARGV[1];

if ( (! $gpept) or (! $species) ) {
    die "Format: gpept.to.prot2gene2sp.pl [GenPept] [species_name] > [protein-to-gene-species.tsv.txt]\n";
}

if ( $species =~ /\W/xms ) { 
    die "Species name \"$species\" has non-word characters\n";
}

open my $GPEPT, '<', $gpept;
while (my $input = <$GPEPT>) {
    chomp $input;
    if ( $input =~ /\A VERSION \s+ (\S+) \s* \z/xms ) {
        $protein = $1;

        # enforce exact, stepwise order of data being read in
        if ($gene) {
            die "Trying to map protein \"$protein\" with incorrect previous gene value \"$gene\"\n";
        }
    }
    elsif ( $input =~ /\A \s+ \/ locus_tag \= \" (\S+) \" \z/xms ) { 
        $gene = $1;

        # again enforce exact, stepwise order of data being read in
        if (! $protein) {
            die "Trying	to map gene \"$gene\" to undefined protein \"$protein\"\n";
        } 

        # if we've made it this far, print a data line
        print "$protein\t$gene\t$species\n";

        # zero the ever-changing $protein and $gene values out
        $protein = q{};
        $gene    = q{};
    }
}
close $GPEPT;

