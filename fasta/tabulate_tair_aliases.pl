#!/usr/bin/env perl

# tabulate_tair_aliases.pl -- Erich Schwarz <ems394@cornell.edu>, 7/21/2015.
# Given ftp://ftp.arabidopsis.org/home/tair/TAIR_Public_Releases/TAIR_Data_20140331/gene_aliases_20140331.txt, make a clean table of TAIR gene IDs to aliases.

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\tAliases";

while (my $input = <>) {
    chomp $input;
    # Sample desired input:
    # AT1G01340       CNGC10  cyclic nucleotide gated channel 10
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M)  G\d+) \t ([^\t]+) \t ([^\t]*) \z/xms ) { 
        my $gene_id = $1;
        my $abbrev  = $2;
        my $name    = $3;
        $name =~ s/\A\s+//;
        $name =~ s/\s+\z//;
        my $full_name = $abbrev;
        if ($name) { 
            $full_name .= " [$name]";
        }
        $data_ref->{'gene_id'}->{$gene_id}->{'full_name'}->{$full_name} = 1;
    }
    elsif ( $input =~ /\A (AT\d+G\d+) \t /xms ) { 
        die "Can't parse input: $input\n";
    }
}

my @gene_ids = sort keys %{ $data_ref->{'gene_id'} };
foreach my $gene_id (@gene_ids) {
    my @full_names = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'full_name'} };
    my $full_name_list = join '; ', @full_names;
    print "$header\n" if $header;
    $header = q{};
    print "$gene_id\t$full_name_list\n";
}

