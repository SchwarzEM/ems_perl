#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $in_data = $ARGV[0];

my $tubul = $ARGV[1];
my $spec  = $ARGV[2];
my $hkeep = $ARGV[3];

my $header = "[Gene, blank]\t[X-axis, blank -- was cell/larvae]\t[Cell_name; was cell TPM]\tTubulin\tSpec\tHkeep";

open my $TUBUL, '<', $tubul;
while ( my $input = <$TUBUL> ) {
    chomp $input;
    $data_ref->{'tubulin'}->{$input} = 1;
}
close $TUBUL;

open my $SPEC, '<', $spec;
while ( my $input = <$SPEC> ) {
    chomp $input;
    $data_ref->{'spec'}->{$input} = 1;
}
close $SPEC;

open my $HKEEP, '<', $hkeep;
while ( my $input = <$HKEEP> ) {
    chomp $input;
    $data_ref->{'hkeep'}->{$input} = 1;
}
close $HKEEP;

open my $IN_DATA, '<', $in_data;
while ( my $input = <$IN_DATA> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) { 
        my $gene   = $1;

        # Note that we are taking *in* three columns -- Gene, Cell_nz in TPM, and Cell_nz/Larvae_nz ratio;
        # we are then swapping the order of Cell_nz and Cell_nz/Larvae_nz so that the latter becomes the x-axis, the former the y-axis;
        # that follows the scheme in Fig. 2B of the LC PNAS paper -- http://www.pnas.org/content/109/40/16246/F2.expansion.html

        my $y_axis = $2;
        my $x_axis = $3;

        # Print this *only* once, at the very start.
        print "$header\n" if $header;
        $header = q{};

        # "[Gene, blank]\t[X-axis, blank -- was cell/larvae]\t[Cell_name; was cell TPM]\tTubulin\tSpec\tHkeep";

        if ( exists $data_ref->{'tubulin'}->{$gene} ) {
            print "$gene\t$x_axis\t\t$y_axis\t\t\n";
        }
        elsif ( exists  $data_ref->{'spec'}->{$gene} ) {
            print "$gene\t$x_axis\t\t\t$y_axis\t\n";
        }
        elsif ( exists $data_ref->{'hkeep'}->{$gene} ) {
            print "$gene\t$x_axis\t\t\t\t$y_axis\n";
        }
        else {
            print "$gene\t$x_axis\t$y_axis\t\t\t\n";
        }
    }
    else {
        die "From input data file $in_data, cannot parse: $input\n";
    }
}
close $IN_DATA;

