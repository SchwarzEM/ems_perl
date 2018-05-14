#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( ( $input !~ /\b NA \b /xms ) and ( $input =~ /\A (WBGene\d+\S+) \t (\S+) \t (\S+) /xms ) ) { 
        my $gene  = $1;
        my $stat1 = $2;
        my $stat2 = $3;
        $stat1 = reformat_stat($stat1);
        $stat2 = reformat_stat($stat2);
        print "$gene\t$stat1\t$stat2\n";
    }
    else { 
        print "$input\n";
    }
}

sub reformat_stat { 
    my $_stat = $_[0];
    if ( $_stat eq '1' ) { 
        return $_stat;
    }
    elsif ( $_stat =~ /\A 0 \. [1-9]+ \d+ \z/xms ) { 
        $_stat = sprintf "%.3f", $_stat;
        return $_stat;
    }
    else {
        $_stat = sprintf "%.3e", $_stat;
        return $_stat;
    }
}

