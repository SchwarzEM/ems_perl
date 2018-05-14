#!/usr/bin/env perl

use strict;
use warnings;

my $consenses = $ARGV[0];
my $pre_ace   = $ARGV[1];

my %mot2cons = ();

open my $CONS, '<', $consenses or die "Can't open consenses table $consenses: $!";
while (my $input = <$CONS>) { 
    chomp $input;
    if ( $input =~ /\A (\d\-\d+) \s+ ([A-Z]+) \s* \z/xms ) { 
        my $mot  = $1;
        my $cons = $2;
        $mot2cons{$mot} = $cons;
    } 
    else { 
        die "From consenses table $consenses, can't parse input: $input\n";
    }
}
close $CONS or die "Can't close filehandle to consenses table $consenses: $!";

open my $PRE_ACE, '<', $pre_ace or die "Can't open precursor .ace file $pre_ace: $!";
while (my $input = <$PRE_ACE>) {
    chomp $input;
    if ( $input =~ /\A (Consensus \s+ \") (\d\-\d+) (\" \s*) \z/xms ) { 
        my $text1 = $1;
        my $mot   = $2;
        my $text2 = $3;
        if (! exists $mot2cons{$mot} ) { 
            die "From precursor .ace file $pre_ace, can't revise: $input\n";
        }
        my $output = $text1 . $mot2cons{$mot} . $text2;
        print "$output\n";
    }
    else { 
        print "$input\n";
    }
}
close $PRE_ACE or die "Can't close filehandle to precursor .ace file $pre_ace: $!"

