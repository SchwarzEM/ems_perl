#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::LineFit;
use Statistics::Descriptive;

my $data_ref;

@{ $data_ref->{'set'}->{'L3_indivs'}->{'cut_index'} }  = (9..13);
@{ $data_ref->{'set'}->{'L4_indivs'}->{'cut_index'} }  = (14..18);
@{ $data_ref->{'set'}->{'nhr_indivs'}->{'cut_index'} } = (19..23);

my @sets = qw( L3_indivs L4_indivs nhr_indivs );

foreach my $set (@sets) { 
    print "\n";

    foreach my $i (0..4) {
        my $query_index = ${ $data_ref->{'set'}->{$set}->{'cut_index'} }[$i];

        my $query_label = `cut -f $query_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | head --lines=1`;
        chomp $query_label;
        my @query_data  = `cut -f $query_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | tail --lines=+2`;
        chomp @query_data;
        $query_label =~ s/_nz\z//; 
        $query_label =~ s/indiv/i/;
        print "$query_label";

        my @j_vals = grep { $_ != $i } (0..4);
        my $data_count = 0;
        foreach my $j (@j_vals) { 
            my $target_index = ${ $data_ref->{'set'}->{$set}->{'cut_index'} }[$j];
            my @target_data = `cut -f $target_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | tail --lines=+2`;
            chomp @target_data;
            @{ $data_ref->{'set'}->{$set}->{'data_index'}->{$j} } = @target_data;
            if ( $data_count > 0 ) { 
                my $new_data_count = @{ $data_ref->{'set'}->{$set}->{'data_index'}->{$j} };
                if ( $new_data_count != $data_count ) {
                    die "Can't parse data because arrays differ in length\n";
                }
            }
            if ( $data_count == 0 ) { 
                $data_count = @{ $data_ref->{'set'}->{$set}->{'data_index'}->{$j} };
            }
        }

        # Switch to zero-based counting of all the elements in the four non-i arrays:
        $data_count--;

        # Build an array whose elements are the mean of the non-i arrays:
        my @mean_of_non_i = ();
        foreach my $datum_number (0..$data_count) {
            my @data_at_number = ();
            foreach my $j (@j_vals) {
                push @data_at_number, ${ $data_ref->{'set'}->{$set}->{'data_index'}->{$j} }[$datum_number];
            }
            my $stat1 = Statistics::Descriptive::Full->new();
            $stat1->add_data(@data_at_number);
            my $data_mean = $stat1->mean();
            push @mean_of_non_i, $data_mean;
        }

        # Having generated it, get r^2 for $i versus the four non-$i:
        my $lineFit = Statistics::LineFit->new();
        $lineFit->setData (\@query_data, \@mean_of_non_i) or die "Invalid data";
        my ($intercept, $slope) = $lineFit->coefficients();
        defined $intercept or die "Can't fit line if x values are all equal";
        my $rSquared = $lineFit->rSquared();
        $rSquared = sprintf("%.2f", $rSquared);
        print "\t$rSquared\n";
    }
}
print "\n";


