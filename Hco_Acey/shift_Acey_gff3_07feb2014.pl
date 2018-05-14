#!/usr/bin/env perl

use strict;
use warnings;

my %scaf2shift = (
    Acey_s0396_scaf => 1,
    Acey_s0536_scaf => 1,
    Acey_s0572_scaf => 1,
    Acey_s0972_scaf => 1510,
    Acey_s1198_scaf => 1673,
    Acey_s1653_scaf => 1,
);

while (my $input = <>) { 
    chomp $input;

    # Sample input:
    # Acey_s0396_scaf	AUGUSTUS	mRNA	9940	10188	.	-	.	ID=Acey_s0396.g664.t1;geneID=Acey_s0396.g664

    if ( $input =~ /\A (\S+) \t ([^\t]* \t [^\t]*) \t (\d+) \t (\d+) \t (.*) \z/xms ) { 
        my $scaf     = $1;
        my $text1    = $2;
        my $start_nt = $3;
        my $stop_nt  = $4;
        my $text2    = $5;
        if (exists $scaf2shift{$scaf} ) { 
            my $shift = $scaf2shift{$scaf};
            $start_nt -= $shift;
            $stop_nt  -= $shift;
        }
        $input = "$scaf\t$text1\t$start_nt\t$stop_nt\t$text2";
    }
    # If a pattern match is not made, the script silently echoes back exactly what it was given by default.
    print "$input\n";
}

