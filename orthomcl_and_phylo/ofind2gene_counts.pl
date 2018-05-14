#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ([A-Z]+\d+) \(\d+ [ ] genes,\d+ [ ] taxa\)[:] \t ([^\t]+) \z/xms ) { 
        my $o_group   = $1;
        my $o_members = $2;

        if ( exists $data_ref->{'o_group'}->{$o_group} ) {
            die "Redundant orthology group (\"$o_group\") in: $input\n";
        }

        while ( $o_members =~ / \S+ \( ( [^\s\(\)]+ ) \) /xmsg ) { 
            my $taxon = $1;
            $data_ref->{'o_group'}->{$o_group}->{'taxon'}->{$taxon}++;
            $data_ref->{'obs_taxa'}->{$taxon} = 1;
        }
    }
    else { 
        die "Cannot parse input line: $input\n";
    }
}

my @obs_taxa       = sort keys %{ $data_ref->{'obs_taxa'} };
my $obs_taxon_list = join "\t", @obs_taxa;
my $header         = "Orthology_group\t$obs_taxon_list";

my @o_groups = sort keys %{ $data_ref->{'o_group'} };

foreach my $o_group (@o_groups) {
    my $output = $o_group;
    foreach my $taxon (@obs_taxa) {
        my $gene_count = 0;
        if ( exists $data_ref->{'o_group'}->{$o_group}->{'taxon'}->{$taxon} ) {
            $gene_count = $data_ref->{'o_group'}->{$o_group}->{'taxon'}->{$taxon} ;
        }
        $output = "$output\t$gene_count";
    }
    print "$header\n" if $header;
    $header = q{};
    print "$output\n";
}


