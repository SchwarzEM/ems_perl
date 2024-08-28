#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $ringers = q{};
my $annots  = q{};

my %seen = ();

$ringers = $ARGV[0] if $ARGV[0];
$annots  = $ARGV[1] if $ARGV[1];

if ( (! $ringers ) or (! $annots ) ) {
    die "Format: get_anomalies_27aug2024.pl [ringers list] [annots]\n";
}

open my $RINGERS, '<', $ringers;
while ( my $input = <$RINGERS> ){
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $seen{$input} = 1;
    }
    else {
        die "From ringers file $ringers, cannot parse: $input\n";
    }
}
close $RINGERS;

open my $ANNOTS, '<', $annots;
while ( my $input = <$ANNOTS> ) {
    chomp $input;
    if ( $input =~ / Parent=(\S+)\.t\d+\t /xms ) {
        my $gene = $1;
        if ( ( exists $seen{$gene} ) and ( $input !~ /Parent=\S+.t\d+ \t gene_id\s\"WBGene\d+\"/xms ) ) {
            print "$input\n";
        }
    }
}
close $ANNOTS;
