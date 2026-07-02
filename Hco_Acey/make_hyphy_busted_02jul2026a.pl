#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $infile = <> ) {
    chomp $infile;
    # Sample input:
    # pangene_1003.dna.mafft_aln.trimal_gpout.fa
    if ( $infile =~ /\A (pangene_\d+\.dna\.mafft_aln\.trimal_gpout)\.fa \z/xms ) {
        my $stem  = $1;
        my $json  = "$stem.busted.json";
        my $log   = "$stem.busted.log.txt";
        my $err   = "$stem.busted.err.txt";

        my $comm1 =   "hyphy CPU=16 busted --alignment $infile --tree orthofinder.1.tree.nwk"
                    . " --error-sink Yes --grid-size 500 --starting-points 5 --kill-zero-lengths No"
                    . " --output $json 1>$log 2>$err ;"
                    ;

        my $comm2 = "if [[ -f errors.log ]]; then";
        my $comm3 = "    cat errors.log >> $err";
        my $comm4 = "    rm errors.log";
        my $comm5 = "fi";

        print "$comm1\n";
        print "$comm2\n";
        print "$comm3\n";
        print "$comm4\n";
        print "$comm5\n";
    }
    else {
        die "Cannot parse infile: $infile\n";
    }
}
