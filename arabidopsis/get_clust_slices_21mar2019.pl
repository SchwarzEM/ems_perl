#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix = q{};
my $infile = q{};

$prefix = $ARGV[0] if $ARGV[0];
$infile = $ARGV[1] if $ARGV[1];

foreach my $i (0..13) {
    my $j = ($i + 1);
    my @inputs = `cat $infile | cut -f $j ;`;
    my $outfile = q{};
    my $OUTFILE;
    foreach my $input (@inputs) {
        chomp $input;

        # Sample input, first line:  C0 (5753 genes)
        if ( $input =~ /\A C\d+ \s+ \( \d+ \s+ genes \) \s* \z/xms) {
            $outfile = $input;

            $outfile =~ s/\s*\z//;
            $outfile =~ s/\s+/_/g;
            $outfile =~ s/\(/_/g;
            $outfile =~ s/\)/_/g;
            $outfile =~ s/[_]+/_/g;
            $outfile =~ s/_\z//;
            $outfile = "$prefix.$outfile.txt";

            open $OUTFILE, '>', $outfile;
        }

        elsif ( $input ne 'Genes' ) {
            print $OUTFILE "$input\n";
        }
    }
    close $OUTFILE;
}
