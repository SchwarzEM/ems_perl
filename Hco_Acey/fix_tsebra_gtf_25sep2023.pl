#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $input  = q{};
my $prefix = q{};

$input  = $ARGV[0] if $ARGV[0];
$prefix = $ARGV[1] if $ARGV[1];

if (! $input ) {
    die "Format: fix_tsebra_gtf_25sep2023.pl [input GFF] [optional gene/tx-name prefix] > [output GFF]\n";
}

if ( $prefix =~ /\A \S+ \z/xms ) {
    $prefix = "$prefix.";
}

open my $INPUT, '<', $input;
while (my $input = <$INPUT>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t/xms ) {
        my $seq = $1;
        $input =~ s/TSEBRA_/$seq.$prefix/g;
        if ( $input =~ /TSEBRA/xms ) {
            die "Failed to fully correct input: $input\n";
        }
        print "$input\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
close $INPUT;
