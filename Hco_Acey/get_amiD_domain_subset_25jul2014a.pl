#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $namelist = $ARGV[0];
my $domains  = $ARGV[1];

my %goodprot = ();

open my $GOOD, '<', $namelist;
while (my $input = <$GOOD>) {
    chomp $input;
    if ( $input =~ /\A \S+ \t (\S+) \z/xms ) { 
        my $prot = $1;
        $goodprot{$prot} = 1;
    }
    else {
        die "From good protein list $namelist, can't parse input line: $input\n";
    }
}
close $GOOD;

open my $DOMS, '<', $domains;
while (my $input = <$DOMS>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (([^\s\/]+) \/ \d+ [-] \d+) /xms ) { 
            my $name = $1;
            my $prot = $2;
            if ( exists $goodprot{$prot} ) { 
                print "$name\n";
            }
        }
        else {
            die "From domain list $domains, can't parse input line: $input\n";
        }
    }
}

