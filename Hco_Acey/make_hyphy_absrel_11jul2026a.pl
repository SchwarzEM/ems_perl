#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $infile = <> ) {
    chomp $infile;
    # Sample input:
    # syntelogs_01/pangene_1003.dna.mafft_aln.trimal_aut1.fa
    if ( $infile =~ /\A ((\S+)\/pangene_\d+\.dna\.mafft_aln\.trimal_aut1)\.fa \z/xms ) {
        my $stem  = $1;
        my $dir   = $2;
        my $json  = "$stem.absrel.json";
        my $log   = "$stem.absrel.log.txt";
        my $err   = "$stem.absrel.err.txt";

        my $comm1 =   "hyphy CPU=16 absrel --alignment $infile --branches Anhui --tree $dir/orthofinder.sp.tree.nwk"
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
