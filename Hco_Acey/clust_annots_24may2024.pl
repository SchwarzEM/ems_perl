#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my @groups = ();
my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+ \t \S+) \t (\S+) \z/xms ) {
        my $group = $1;
        my $rep   = $2;
        push @groups, $group;
        push @{ $data_ref->{'group'}->{$group}->{'reps'} }, $rep;
    }
}

@groups = uniq(@groups);

foreach my $group (@groups) {
    my @reps = @{ $data_ref->{'group'}->{$group}->{'reps'} };
    my $rep_text = join ', ', @reps;
    print "$group\t$rep_text\n";
}

