#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $treefam_id = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (TF\d+) \s* \z/xms ) { 
        $treefam_id = $1;
    }
    elsif ( $input =~ /\A [#] END \s* \z/xms ) {
        $treefam_id = q{};
    }
    elsif ( $treefam_id and ( $input =~ /\A > (\S+) /xms ) ) {
        my $prot_id = $1;
        $data_ref->{'treefam_id'}->{$treefam_id}->{'prot_id'}->{$prot_id} = 1;
    }
}

my @treefam_ids = sort keys %{ $data_ref->{'treefam_id'} };
foreach my $treefam_id1 (@treefam_ids) {
    my @prot_ids = sort keys %{ $data_ref->{'treefam_id'}->{$treefam_id1}->{'prot_id'} };
    foreach my $prot_id (@prot_ids) {
        print "$treefam_id1\t$prot_id\n";
    }
}


