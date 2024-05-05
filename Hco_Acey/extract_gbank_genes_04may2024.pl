#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $read    = 0;
my $gene_id = q{};
my $gb_id   = q{};

while ( my $input = <> ) {
    # Sample input:
    #    gene            23529..25090
    # [or]
    #    gene            complement(24619..28447)
    #                    /locus_tag="QR680_007293"
    #                    /note="Sherm_I.g354"
    chomp $input;

    if ( ( $read ) and ( $gb_id ) ) {
        if ( $gene_id ) {
            die "Redundant native gene IDs: $gene_id [and] $input\n";
        }
        elsif ( $input =~ /\A \s+ \/note=\"(\S+)\"/xms  ) {
            $gene_id = $1;
        }
        else {
            die "Cannot parse what should be a native gene ID: $input\n";
        }

        print "$gene_id\t$gb_id\n";

        $read    = 0;
        $gene_id = q{};
        $gb_id   = q{};
    }
    elsif ( ( $read ) and (! $gene_id ) ) {
        if ( $input =~ /\A \s+ \/locus_tag=\"(\S+)\" /xms ) {
            $gb_id   = $1;
            $gene_id = q{};
        }
        else {
            die "Cannot parse what should be a GenBank gene ID: $input\n";
        }
    }
    elsif ( (! $read ) and ( $input =~ /\A \s+ gene \s+ /xms ) ) {
        $read    = 1;
        $gene_id = q{};
        $gb_id   = q{};
    }
}
