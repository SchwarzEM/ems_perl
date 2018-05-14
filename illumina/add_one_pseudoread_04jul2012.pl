#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A WBGene\d+\S+ \t (?: \d+ \t )+ \d+ \z/xms ) { 
        if ( $input =~ /\A \t/xms ) { 
            print "$input\n";
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
    elsif ( $input =~ /\A (WBGene\d+\S+) \t ((?: \d+ \t )+ \d+) \z/xms ) { 
        my $wbgene            = $1;
        my $readcount_inline  = $2;
        my @readcounts_input  = split /\t/, $readcount_inline;
        my @readcounts_output = ();

        foreach my $readcount (@readcounts_input) { 
            $readcount++;
            push @readcounts_output, $readcount;
        }

        my $readcount_outline = join "\t", @readcounts_output;

        # "\t" not '\t' !
        my $output = $wbgene . "\t" . $readcount_outline;

        print "$output\n";
    }
    else {
        die "*Really* can't parse input: $input\n";
    }
}

