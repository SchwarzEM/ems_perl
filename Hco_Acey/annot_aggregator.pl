#!/usr/bin/env perl

# annot_aggregator.pl -- Erich Schwarz <ems394@cornell.edu>, 9/30/2014.
# Purpose: given, as an argument, a user-specified data-column header title ("Title"), and a simple table with "gene name \t some annotation", generate an aggregated table with the header "Gene\tTitle", followed by lines where each gene is followed by *all* its annotations.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $data_title  = q{};
my $data_table = q{};

my $help;

GetOptions ( 'title=s' => \$data_title,
             'data=s'  => \$data_table,
             'help'    => \$help, 
);

if ( $help or (! $data_title) or (! $data_table) ) { 
    die "Format: annot_aggregator.pl\n",
        "    --title|-t  [user-specified data-column header title]\n",
        '    --data|-d   [simple table with "gene name \t some annotation"]', "\n",
        "    --help|-h   [print this message]\n",
        ;
}

my $header = "Gene\t$data_title";

open my $DATA, '<', $data_table;
while (my $input = <$DATA>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) { 
        my $gene  = $1;
        my $annot = $2;
        $data_ref->{'gene'}->{$gene}->{'annot'}->{$annot} = 1;
    }
    else {
        die "From data table $data_table, can't parse input: $input\n";
    }
}
close $DATA;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    print "$header\n" if $header;
    $header = q{};
    my @annots = sort keys %{ $data_ref->{'gene'}->{$gene}->{'annot'} };
    my $annot_text = join '; ', @annots;
    print "$gene\t$annot_text\n";
}


