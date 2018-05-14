#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %pref2spec = (
    RumMG_ => 'Sheep rumen metagenome',
    StMG1_ => 'Human stool metagenome (stool_uni_mh.fa)',
    StMG2_ => 'Human stool metagenome (stool_nr.fa)',
);

my @proteomes = @ARGV;
my $sec_amids = pop @proteomes; 
my %seen      = ();

open my $SEC, '<', $sec_amids;
while (my $input = <$SEC>) {
    chomp $input;
    if ( $input =~ /\A > Metg_ ([^\s\/]+) \/ \d+ [-] \d+ /xms ) {
        my $protein = $1;
        $seen{$protein} = 1;
    }
}
close $SEC;

foreach my $proteome (@proteomes) {
    open my $PROT, '<', $proteome;
    while (my $input = <$PROT>) { 
        chomp $input;
        if ( $input =~ /\A > ( ([A-Z] [A-Za-z0-9]+ [_]) (\S+) ) /xms ) { 
            my $part_name   = $1;
            my $metg_prefix = $2;
            my $protein     = $3;
            my $full_name   = 'Metg_' . $part_name;
            if ( ( exists $seen{$part_name} ) and ( exists $pref2spec{$metg_prefix} ) ) {
                my $species = $pref2spec{$metg_prefix};
                print "$full_name\tMetagenomes\t$species\t$protein\t$protein\n";
            }
        }
    }
    close $PROT;
}

