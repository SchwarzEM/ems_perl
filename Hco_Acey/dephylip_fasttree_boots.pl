#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix  = $ARGV[0];
my $synos   = $ARGV[1];
my $inboots = $ARGV[2];

my %syn2real = ();

if ( (! $prefix) or (! $synos) or (! $inboots) ) { 
    die "Format: dephylip_fasttree_boots.pl [phylip synonym prefix] [synonym table] [input bootstrap table file] > [output, renamed bootstrap table file]\n";
}

open my $SYNOS, '<', $synos;
while (my $input = <$SYNOS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $phylip_name = $1;
        my $real_name   = $2;
        if ( exists $syn2real{$phylip_name} ) {
            die "Redundant use of PHYLIP name: $phylip_name\n";
        }
        $syn2real{$phylip_name} = $real_name;
    }
    else {
        die "From synonym table file $synos, can't parse input line: $input\n";
    }
}
close $SYNOS;

open my $INBOOTS, '<', $inboots;
while (my $input = <$INBOOTS>) {
    chomp $input;
    while ( $input =~ /\A (.*?) ($prefix\w+) (.*) \z/xmsg ) {
        my $front_text  = $1; 
        my $phylip_name = $2;
        my $back_text   = $3;
        if (! exists $syn2real{$phylip_name} ) {
            die "In bootstrap tree file $inboots, failed to parse putative PHYLIP synonym $phylip_name\nContext line: $input\n";
        }
        my $real_name = $syn2real{$phylip_name};
        $input = $front_text . $real_name . $back_text;
    }
    print "$input\n";
}
close $INBOOTS;




