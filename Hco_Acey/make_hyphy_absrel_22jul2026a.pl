#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $in_list = q{};
my $taxon   = q{};

$in_list = $ARGV[0] if $ARGV[0];
$taxon   = $ARGV[1] if $ARGV[1];

if ( (! $in_list ) or (! $taxon ) ) {
    die "Format: make_hyphy_absrel_22jul2026a.pl [input file list] [target taxon] > [HyPhy aBSREL branch-specific commands]\n";
}

open my $IN_LIST, '<', $in_list;
while ( my $infile = <$IN_LIST> ) {
    chomp $infile;
    # Sample input:
    # syntelogs_01/pangene_1003.dna.mafft_aln.trimal_aut1.fa
    if ( $infile =~ /\A ((\S+)\/pangene_\d+\.dna\.mafft_aln\.trimal_aut1)\.fa \z/xms ) {
        my $stem  = $1;
        my $dir   = $2;
        my $json  = "$stem.$taxon.absrel.json";
        my $log   = "$stem.$taxon.absrel.log.txt";
        my $err   = "$stem.$taxon.absrel.err.txt";

        my $comm1 =   "hyphy CPU=16 absrel --alignment $infile --branches $taxon --tree $dir/orthofinder.sp.tree.nwk"
                    . " --error-sink Yes --grid-size 500 --starting-points 5 --kill-zero-lengths No"
                    . " --output $json 1>$log 2>$err ;"
                    ;

        my $comm2 = "if [[ -f errors.log ]]; then";
        my $comm3 = "    cat errors.log >> $err";
        my $comm4 = "    echo           >> $err";
        my $comm5 = "    rm errors.log";
        my $comm6 = "fi";

        print "$comm1\n";
        print "$comm2\n";
        print "$comm3\n";
        print "$comm4\n";
        print "$comm5\n";
        print "$comm6\n";
    }
    else {
        die "Cannot parse infile: $infile\n";
    }
}
close $IN_LIST;
