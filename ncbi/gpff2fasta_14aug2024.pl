#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $protein = q{};
my $gene    = q{};
my $read    = 0;
my $seq     = q{};

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
        print ">$protein\tgene=$gene\n";

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
