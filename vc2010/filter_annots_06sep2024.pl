#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $rejects = q{};
my $annots  = q{};

$rejects = $ARGV[0] if $ARGV[0];
$annots  = $ARGV[1] if $ARGV[1];

my %reject  = ();

if ( (! $rejects ) or (! $annots ) ) {
    die "Format: filter_annots_06sep2024.pl [rejects] [annots] > [filtered annots]\n";
}

open my $REJECTS, '<', $rejects;
while ( my $input = <$REJECTS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \b/xms ) {
        my $gene = $1;
        $reject{$gene} = 1;
    }
    else {
        die "From rejects $rejects, cannot parse: $input\n";
    }
}
close $REJECTS;

open my $ANNOTS, '<', $annots;
while ( my $input = <$ANNOTS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t /xms ) {
        my $gene = $1;
        if ( (! exists $reject{$gene} ) and ( $gene ne 'Gene' ) ) {
            print "$input\n";
        }
    }
    else {
        die "From annots $annots, cannot parse: $input\n";
    }
}
close $ANNOTS;
