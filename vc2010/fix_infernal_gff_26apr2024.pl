#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ( \S+ \t cmscan \t (\S+) (?:\t[^\t]*){5} ) \t (\S+) \z/xms ) {
        my $annot1   = $1;
        my $rfam_tag = $2;
        my $annot2   = $3;

        my $e_value  = q{};
        my $rfam_acc = q{};
        my $desc     = q{};
        my $annot3   = q{};
        
        if ( $annot2 =~ /\A evalue = ([^;\s]+) ; \S+ ; mdlaccn = ([^;\s]+) ; \S+ ; desc=\d+_\d+_(\S+) \z/xms ) { 
            $e_value  = $1;
            $rfam_acc = $2;
            $desc     = $3;
            $annot3 = "evalue=$e_value;mdlaccn=$rfam_acc|$rfam_tag;desc=$desc";
        }
        else {
            die "Cannot parse annot2: $annot2\n";
        }

        print "$annot1\t$annot3\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
