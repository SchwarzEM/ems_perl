#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $washu_cds2gene = q{};
my $washu2umass    = q{};

$washu_cds2gene = $ARGV[0] if $ARGV[0];
$washu2umass    = $ARGV[1] if $ARGV[1];

my %cds2gene = ();

open my $WASHU_CDS2GENE, '<', $washu_cds2gene;
while (my $input = <$WASHU_CDS2GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $cds  = $1;
        my $gene = $2;
        $cds2gene{$cds} = $gene;
    }
    else {
        die "In $washu_cds2gene, cannot parse: $input\n";
    }
}
close $WASHU_CDS2GENE;

open my $WASHU2UMASS, '<', $washu2umass;
while (my $input = <$WASHU2UMASS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $cds   = $1;
        my $gene2 = $2;
        if (! exists $cds2gene{$cds} ) {
            die "In $washu2umass, cannot identify a gene for cds: $cds\n";
        }
        my $gene1 = $cds2gene{$cds};
        print "$gene1\t$gene2\n";
    }
    else {
        die "In $washu2umass, cannot parse: $input\n";
    }
}
close $WASHU2UMASS;

