#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $wormpep = q{};
my $cds_dna = q{};

$wormpep = $ARGV[0] if $ARGV[0];
$cds_dna = $ARGV[1] if $ARGV[1];

if ( (! $wormpep ) or (! $cds_dna ) ) {
    die "Format: wormpep2cds_headers.pl [wormpep] [cds_dna] > [cds_dna_with_wormpep_headers] ;\n";
}

my %cds2header = ();

open my $WP, '<', $wormpep;
while (my $input = <$WP> ) {
    chomp $input;
    if ( $input =~ /\A [>] ((\S+) .*) \z/xms ) {
        my $header = $1;
        my $cds    = $2;
        $cds2header{$cds} = $header;
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "In wormpep $wormpep, cannot parse header: $input\n";
    }
}
close $WP;

open my $CDS, '<', $cds_dna;
while (my $input = <$CDS> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) /xms ) {
        my $cds    = $1;
        my $header = q{};

        if ( exists $cds2header{$cds} ) {
            $header = $cds2header{$cds};
            delete $cds2header{$cds};
        }
        else {
            die "In CDS file $cds_dna, cannot map to wormpep header: $input\n";
        }

        print ">$header\n";
    }
    elsif ( $input =~ /\A [>] /xms ) {
       	die "In	CDS file $cds_dna, cannot parse header: $input\n";
    }
    else {
        print "$input\n";
    }
}
close $CDS;

