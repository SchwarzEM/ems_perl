#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile  = q{};
my $phobius = q{};

$infile  = $ARGV[0] if $ARGV[0];
$phobius = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $phobius ) ) {
    die "Format: bool_phobius_18oct2025.pl [gene-Phobius annots.] [spec. Phobius annot.] > [Boolean table for spec. Phobius annot.]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    if ( $input =~ /\A Gene \t/xms ) {
        $input = "Gene\t$phobius";
    }
    elsif ( ( $input !~ /\A Gene \t/xms ) and ( $input =~ /\A (\S+) \t/xms ) ) {
        my $gene = $1;
        if ( ( $input =~ /\A \S+ \t .* \b $phobius [;] .* /xms ) or ( $input =~ /\A \S+ \t .* \b $phobius \z/xms ) ) {
            $input = "$gene\ttrue";
        }
        else {
            $input = "$gene\tfalse";
        }
    }
    else {
        die "From input file $infile, cannot parse: $input\n"
    }
    print "$input\n";
}
close $INFILE;


