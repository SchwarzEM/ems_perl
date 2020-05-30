#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $elepep2wormpep = q{};
my $wormpep2repeat = q{};

$wormpep2repeat = $ARGV[0] if $ARGV[0];
$elepep2wormpep = $ARGV[1] if $ARGV[1];

my $data_ref;

if ( (! $elepep2wormpep) or (! $wormpep2repeat ) ) {
    die "Format: wormpep2rep_30may2020.pl [wormpep2repeat] [elepep2wormpep] > [elements_to_disqualify] ;\n";
}

open my $WORMPEP2REPEAT, '<', $wormpep2repeat ;
while (my $input = <$WORMPEP2REPEAT>) {
    chomp $input;
    if ( $input =~ /\A \S+ \s+ \S+ \s+ (\S+) \s+ \S.+ \z/xms ) {
        my $wormpep = $1;
        $data_ref->{'rep_wpep'}->{$wormpep} = 1;
    }
    else {
        die "From wormpep2repeat $wormpep2repeat, cannot parse: $input\n";
    }
}
close $WORMPEP2REPEAT;

open my $ELEPEP2WORMPEP, '<', $elepep2wormpep;
while (my $input = <$ELEPEP2WORMPEP>) {
    chomp $input;
    if ( $input =~ /\A (\S+) _\d+ \t (\S+) \z/xms ) {
        my $element = $1;
        my $wormpep = $2;
        if (! exists $data_ref->{'rep_wpep'}->{$wormpep} ) {
            $data_ref->{'element'}->{$element} = 1;
        }
    }
    else {
        die "From elepep2wormpep $elepep2wormpep, cannot parse: $input\n";
    }
}
close $ELEPEP2WORMPEP;

my @elements = sort keys %{ $data_ref->{'element'} };
foreach	my $element (@elements) {
    print "$element\n";
}

