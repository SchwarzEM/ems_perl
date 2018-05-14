#!/usr/bin/env perl

use strict;
use warnings;

my $phobius_file = $ARGV[0];
my $flp_mot_file = $ARGV[1];

my %has_sigp = ();

open my $PHOBIUS, '<', $phobius_file or die "Can't open phobius file $phobius_file: $!";
while (my $input = <$PHOBIUS>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ SigP \z/xms ) { 
        my $cds = $1;
        $has_sigp{$cds} = 1;
    }
}
close $PHOBIUS or die "Can't close handle to phobius file $phobius_file: $!";

open my $FLP, '<', $flp_mot_file or die "Can't open flp motif file $flp_mot_file: $!";
while (my $input = <$FLP>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ .+ \z/xms ) {
        my $cds = $1;
        if ( exists $has_sigp{$cds} ) { 
            print "$input\n";
        }
    }
}
close $FLP or die "Can't close handle to flp motif file $flp_mot_file: $!";

