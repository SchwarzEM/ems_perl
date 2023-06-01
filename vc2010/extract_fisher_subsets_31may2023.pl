#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $list = q{};
$list    = $ARGV[0] if $ARGV[0];

if (! $list ) {
    die "Format: extract_fisher_subsets_31may2023.pl [list of BEDtools fisher output files] > [TSV: motif name / ratio / two-tailed p-value ]\n";
}

open my $LIST, '<', $list;
while (my $infile = <$LIST>) {
    chomp $infile ;

    my $stem;
    my $basename = basename($infile);

    if ( $basename =~ /\A (\S+) \.fisher\.\S+\.txt \z/xms ) {
        $stem = $1;
    }
    else {
        die "Cannot parse file: $infile\n";
    }

    my $data = `cat $infile | tail --lines=1 | cut -f 3-4 ;`;
    chomp $data;
    my $p_value = q{};
    my $ratio   = q{};
    if ( $data =~ /\A (\S+) \t (\S+) \z/xms ) {
        $p_value = $1;
        $ratio   = $2;
    }
    else {
        die "Cannot parse bedtools fisher data: $data\n";
    }
    print "$stem\t$ratio\t$p_value\n";
}
close $LIST;
