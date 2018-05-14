#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $asp_file  = $ARGV[0];
my $aspr_file = $ARGV[1];

if ( (! $asp_file) or (! $aspr_file) ) { 
    die "Format: generate_nam_asp_or_aspr_table_25mar2014.pl [Nam ASP gene list] [Nam ASPR gene list]\n";
}

my $header = "Synonym\tGene\tType";

my %synonyms = ( 
    'NECAME_09334' => 'Nam-ASP-1',
    'NECAME_11281' => 'Nam-ASP-2',
    'NECAME_15644' => 'Nam-ASP-6',
    'NECAME_01333' => 'Nam-ASP-7',
);

my @output_lines = ();

open my $ASP, '<', $asp_file;
while (my $input = <$ASP>) { 
    chomp $input;
    if ( $input =~ /\A (NECAME_(\d{5})) \z/xms ) { 
        my $cds   = $1;
        my $label = $2;
        my $gene  = 'Nam-ASP-g' . $label;
        if ( exists $synonyms{$cds} ) {
            $gene = $synonyms{$cds};
        }
        my $output_line = "$gene\t$cds\tASP";
        push @output_lines, $output_line;
    }
    else {
        die "From asp_file $asp_file, can't parse: $input\n";
    }
}
close $ASP;

open my $ASPR, '<', $aspr_file;
while (my $input = <$ASPR>) {
    chomp $input;
    if ( $input =~ /\A (NECAME_(\d{5})) \z/xms ) {
        my $cds   = $1;
        my $label = $2;
        my $gene  = 'Nam-ASPR-g' . $label;
        my $output_line = "$gene\t$cds\tASPR";
        push @output_lines, $output_line;
    }
    else {
        die "From asp_file $asp_file, can't parse: $input\n";
    }
}
close $ASPR;

@output_lines = sort @output_lines;

if (@output_lines) {
    print "$header\n";
    foreach my $output_line (@output_lines) {
        print "$output_line\n";
    }
}

