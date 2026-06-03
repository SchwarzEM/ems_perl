#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $protein  = q{};
my $alt_prot = q{};
my $gene     = q{};
my $read     = 0;
my $seq      = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A VERSION \s+ (\S+) /xms ) {
        $protein = $1;
        $gene    = q{};
        $read    = 0;
        $seq     = q{};
    }
    elsif ( $input =~ / gene [=] ["] (\S+) ["] /xms ) {
        $gene = $1;
        $read = 0;
        $seq  = q{};
    }
    # Not all GenPept files have 'gene' tags; others have 'locus_tag'.
    elsif ( $input =~ / locus_tag [=] ["] (\S+) ["] /xms ) {
        $gene = $1;
        $read = 0;
        $seq  = q{};
    }
    # Alternative protein names are sometimes given in annotations like
    #     'ID:CAUJ4_0100000100.t1.cds', 'ID:CBOVI.g1.t1.CDS2', 'ID:CAEBRE_00001.1-cds2', or 'ID:Sp36v6_10000100.t1.cds2'.

    elsif ( $input =~ / ID [:] (\S+?) [;] /xms ) {
        $alt_prot = $1;
        $read     = 0;
        $seq      = q{};

        # Trim off various suffixes:
        if ( $alt_prot =~ /\A\S+(?:-|\.)cds[\d]*\z/ ) {
            $alt_prot =~ s/\A(\S+)(?:-|\.)cds[\d]*\z/$1/;
        }
        elsif ( $alt_prot =~ /\A\S+(?:-|\.)CDS[\d]*\z/   ) {
            $alt_prot =~ s/\A(\S+)(?:-|\.)CDS[\d]*\z/$1/;
        }
    }

    elsif ( $input =~ / ORIGIN /xms ) {
        $read = 1;
    }
    elsif ( $read and ( $input !~ /\A \/\/ /xms ) and ( $input =~ / \S /xms ) ) {
        $input =~ s/\d//g;
        $input =~ s/\s//g;
        $input =~ tr/[a-z]/[A-Z]/;
        $seq   = $seq . $input;
    }
    elsif ( $input =~ /\A \/\/ /xms ) {
        $read = 0;
        print ">$protein\t$alt_prot\tgene=$gene\n";

        my @output_lines = unpack("a60" x (length($seq)/60 + 1), $seq);
        foreach my $output_line (@output_lines) { 
            if ($output_line =~ /\S/) { 
                print "$output_line\n";
            }
        }

        $protein = q{};
        $gene    = q{};
        $seq     = q{};
    }
}
