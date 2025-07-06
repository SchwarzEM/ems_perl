#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gene_bed   = q{};
my $pixy_stats = q{};

$gene_bed   = $ARGV[0] if $ARGV[0];
$pixy_stats = $ARGV[1] if $ARGV[1];

if ( (! $gene_bed ) or (! $pixy_stats ) ) {
    die "Format: annot_pixy_05jul2025a.pl [gene BED] [pixy stats]\n";
}

my $data_ref;

my $header = "Gene\tdxy_1\tdxy_2\tRatio";

open my $GENES, '<', $gene_bed;
while ( my $input = <$GENES> ) {
    chomp $input;
    # Sample input:
    # Necator_chrI    13891   14485   .       .       -       AUGUSTUS        gene    .       ID=Necator_chrI.g137
    if ( $input =~ /\A (\S+ \t \d+ \t \d+) \t .* \t ID=(\S+) \s* \z/xms ) {
        my $coords = $1;
        my $gene   = $2;
        $data_ref->{'coords'}->{$coords}->{'gene'}->{$gene} = 1;
    }
    else {
        die "From gene BED $gene_bed, cannot parse: $input\n";
    }
}
close $GENES;

open my $PIXY, '<', $pixy_stats;
while ( my $input = <$PIXY> ) {
    chomp $input;
    my $coords = q{};
    my $data   = q{};
    my @genes  = ();
    # Sample input:
    # Necator_chrI    8404666 8408033 0.0     0.2857142857142857      285.714285714286
    if ( $input =~ /\A (\S+ \t \d+ \t \d+) \t (\S+ \t \S+ \t \S+) \s* \z/xms ) {
        $coords = $1;
        $data   = $2;
    }
    else {
        die "From data stats $pixy_stats, cannot parse: $input\n";
    }
    if ( exists $data_ref->{'coords'}->{$coords}->{'gene'} ) {
        @genes = sort keys %{ $data_ref->{'coords'}->{$coords}->{'gene'} };
        foreach my $gene (@genes) {
            print "$header\n" if $header;
            $header = q{};
            print "$gene\t$data\n";
        }
    }
}
close $PIXY;
