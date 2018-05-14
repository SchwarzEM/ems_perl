#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $input_data   = $ARGV[0];
my $gpcrs        = $ARGV[1];
my $housekeeping = $ARGV[2];

if ( (! $ARGV[0]) or (! $ARGV[1]) or (! $ARGV[2]) ) {
    die "Format: awc_tablemaker_25may2015_v01.pl [input data] [GPCR gene list] [housekeeping gene list]\n";
}

my %is_gpcr  = ();
my %is_house = ();

# Replace this: Gene	AWC/larvae	AWC_nz
my $header = "Gene\t\tAWC\tGPCR\tHousekeeping\n";

my $output = q{};

open my $GPCRS, '<', $gpcrs;
while ( my $input = <$GPCRS> ) {
    chomp $input;
    $is_gpcr{$input} = 1;
}
close $GPCRS;

open my $HOUSE, '<', $housekeeping;
while ( my $input = <$HOUSE> ) {
    chomp $input;
    $is_house{$input} = 1;
}
close $HOUSE;

open my $DATA, '<', $input_data;
while ( my $input = <$DATA> ) {
    # Sample input line:
    # WBGene00000001|Y110A7A.10|aap-1        74.20   134.30
    chomp $input;
    print $header if $header;
    $header = q{};
    if ( $input =~ /\A ((WBGene\S+) \t \S+) \t (\S+) \z/xms ) { 
        my $front_line = $1;
        my $wbgene     = $2;
        my $awc_tpm    = $3;
        if ( exists $is_house{$wbgene} ) { 
            $output = $front_line . "\t\t\t$awc_tpm";
        }
        elsif ( $is_gpcr{$wbgene} ) {
            $output = $front_line . "\t\t$awc_tpm\t";
        }
        else {
            $output = $front_line . "\t$awc_tpm\t\t";
        }
        print "$output\n";
    }
    else {
        warn "From input data file $input_data, could not parse: $input\n";
    }
}
close $DATA;


