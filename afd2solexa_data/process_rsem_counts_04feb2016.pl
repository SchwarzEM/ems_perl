#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $readset      = q{};
my $total_reads  = q{};
my $aligned_zero = q{};
my $aligned_one  = q{};
my $aligned_two  = q{};
my $align_rate   = q{};
my $header       = "Data\tTotal_reads\tAligned_0_times\tAligned_exactly_1_time\tAligned_2+_times\tOverall_alignment_rate";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A Input: \s+ (\S+) /xms ) { 
        $readset      = $1;
        $total_reads  = q{};
        $aligned_zero = q{};
        $aligned_one  = q{};
        $aligned_two  = q{};
        $align_rate   = q{};
    }
    elsif ( $readset and ( $input =~ /\A (\d+) \s reads; \s of \s these: \s* \z/xms  ) ) {
        $total_reads = $1;
    }
    elsif ( $readset and ( $input =~ /\A \s+ (\d+ \s \(\S+\)) \s aligned \s 0 \s times \s* \z/xms ) ) { 
        $aligned_zero = $1;
    }
    elsif ( $readset and ( $input =~ /\A \s+ (\d+ \s \(\S+\)) \s aligned \s exactly \s 1 \s time \s* \z/xms ) ) {
        $aligned_one = $1;
    }
    elsif ( $readset and ( $input =~ /\A \s+ (\d+ \s \(\S+\)) \s aligned \s [>] 1 \s times \s* \z/xms ) ) {
        $aligned_two = $1;
    }
    elsif ( $readset and ( $input =~ /\A (\S+) \s overall \s alignment \s rate \s* \z/xms ) ) {
        $align_rate = $1;

        print "$header\n" if $header;
        $header = q{};

        print "$readset\t$total_reads\t$aligned_zero\t$aligned_one\t$aligned_two\t$align_rate\n";

        $readset      = q{};
        $total_reads  = q{};
        $aligned_zero = q{};
        $aligned_one  = q{};
        $aligned_two  = q{};
        $align_rate   = q{};
    }
}

