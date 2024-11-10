#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $list = q{};
my $gff  = q{};

my %listed = ();
my %seen   = ();

$list = $ARGV[0] if $ARGV[0];
$gff  = $ARGV[1] if $ARGV[1];

if ( (! $list ) or (! $gff ) ) {
    die "Format: gff_gene_subset_09nov2024a.pl [gene list] [gff to filter] > [gff with only listed genes]\n";
}

open my $LIST, '<', $list;
while ( my $gene = <$LIST> ) {
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        $listed{$gene} = 1;
    }
    else {
        die "From gene list $list, cannot parse: $gene\n";
    }
}
close $LIST;

open my $GFF, '<', $gff;
while ( my $input = <$GFF> ) {
    chomp $input;
    if ( ( $input =~ / ID[=](\S+)[;] /xms ) and ( exists $listed{$1} ) ) {
        my $gene = $1;
        $seen{$gene} = 1;
        print "$input\n";
    }
    elsif ( $input !~ / ID[=](\S+)[;] /xms ) {
        warn "From GFF $gff, could not parse: $input\n";
    }
}

my @listeds = sort keys %listed;
foreach my $gene (@listeds) {
    if (! exists $seen{$gene} ) {
        warn "Listed but never seen in the GFF: $gene\n";
    }
}
