#!/usr/bin/env perl

use strict;
use warnings;

my $infile = $ARGV[0];
my $wanted = $ARGV[1];

my $gene  = q{};
my $group = q{};

# Just asking for (! $wanted) fails in cases where the argument is, in fact, '0' == 'no'...
# So instead ask for a $wanted that has at least one non-space character.

if (! ( $wanted =~ /\S/xms ) ) { 
    die "Format: kgg2group_08dec2010.pl [input KGG file] [desired group]\n";
}

open my $INFILE, '<', $infile or die "Can't open input KGG file $infile: $!";
while (my $input = <$INFILE>) { 
    chomp $input;
    if ( $input=~ /\A (WBGene\d+) \S* \t (\S+) /xms ) { 
        $gene = $1;
        $group = $2;
        if ( $group eq $wanted ) { 
            print "$gene\n";
        }
    }
}

