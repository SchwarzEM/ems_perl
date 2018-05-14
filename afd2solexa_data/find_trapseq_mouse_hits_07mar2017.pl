#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $mouse_hits = $ARGV[0];
my $tair_hits  = $ARGV[1];

my $cDNA        = q{};
my $mouse_hit   = q{};
my $mouse_e_val = q{};
my $tair_e_val  = q{};

my $data_ref;

open my $MOUSE, '<', $mouse_hits;
while (my $input = <$MOUSE>) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (\S+) \t (\S+) \t .+ \t (\S+) \t \S+ \z/xms ) ) {
        $cDNA        = $1;
        $mouse_hit   = $2;
        $mouse_e_val = $3;
        $data_ref->{'cDNA'}->{$cDNA}->{'mouse_hit'}   = $mouse_hit;
        $data_ref->{'cDNA'}->{$cDNA}->{'mouse_e_val'} = $mouse_e_val;
    }
}
close $MOUSE;

open my $TAIR, '<', $tair_hits;
while (my $input = <$TAIR>) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (\S+) \t \S+ \t .+ \t (\S+) \t \S+ \z/xms ) ) {
        $cDNA       = $1;
        $tair_e_val = $2;
        $data_ref->{'cDNA'}->{$cDNA}->{'tair_e_val'} = $tair_e_val;
    }
}
close $TAIR;

my @cDNAs = sort grep { exists $data_ref->{'cDNA'}->{$_}->{'mouse_e_val'} } keys %{ $data_ref->{'cDNA'} };

foreach my $cDNA (@cDNAs) {
    if ( (! exists $data_ref->{'cDNA'}->{$cDNA}->{'mouse_hit'} ) or (! exists $data_ref->{'cDNA'}->{$cDNA}->{'mouse_e_val'} ) ) {
        die "Failed to map cDNA $cDNA to mouse hit and mouse E-value\n";
    }
    $mouse_hit   = $data_ref->{'cDNA'}->{$cDNA}->{'mouse_hit'};
    $mouse_e_val = $data_ref->{'cDNA'}->{$cDNA}->{'mouse_e_val'};

    $tair_e_val  = 1;
    if ( exists $data_ref->{'cDNA'}->{$cDNA}->{'tair_e_val'} ) {
        $tair_e_val  = $data_ref->{'cDNA'}->{$cDNA}->{'tair_e_val'};
    }

    if ( $mouse_e_val < $tair_e_val ) {
        print "$cDNA\t$mouse_hit\t$mouse_e_val\n";
    }
}

