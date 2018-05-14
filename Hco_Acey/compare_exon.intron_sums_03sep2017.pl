#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Math::GSL::Fit qw /:all/;
use Math::GSL::Statistics qw /:all/;;

my $data_ref;

my $header = "Comparison\tr^2\tSlope\ty-intersect";

my @data_types = ();

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A Orthogroup \t /xms ) { 
        @data_types = split /\t/, $input;
    }
    elsif (@data_types) {
        my @data_values = split /\t/, $input;
        my $entry_no = @data_values;
        $entry_no--;
        foreach my $i (0..$entry_no) {
            push @{ $data_ref->{$data_types[$i]} }, $data_values[$i];
        }
    }
    else {
        die "Cannot parse: $input\n";
    }
}
close $INFILE;

my @key_exon_data_types = ();
@key_exon_data_types    = grep { /exon/ } grep { /_sum\z/ } @data_types;

my @key_intron_data_types = ();
@key_intron_data_types    =  grep { /intron/ } grep { /_sum\z/ } @data_types;

my @key_data_type_refs = (\@key_exon_data_types, \@key_intron_data_types);

# Compare each key data type to all others -- but only once.
foreach my $key_data_type_ref (@key_data_type_refs) {
    my @key_data_types = @{ $key_data_type_ref };
    while (@key_data_types) {
        my $key_data1 = shift @key_data_types;

        foreach my $key_data2 (@key_data_types) {
            my $key_data1_values_ref = $data_ref->{$key_data1};
            my $key_data2_values_ref = $data_ref->{$key_data2};

            my $n = @{$key_data1_values_ref};
            my ($status, $y_intersect, $slope, $cov00, $cov01, $cov11, $ss_resid) 
                = gsl_fit_linear($key_data1_values_ref, 1, $key_data2_values_ref, 1, $n);

            my $pearson = gsl_stats_correlation($key_data1_values_ref, 1, $key_data2_values_ref, 1, $n); 
            my $r2 = ($pearson**2);

            print "$header\n" if $header;
            $header = q{};

            print "$key_data1.vs.$key_data2\t$r2\t$slope\t$y_intersect\n";
        }
    }
}

