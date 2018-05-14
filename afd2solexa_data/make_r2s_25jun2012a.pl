#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::LineFit;

my $data_ref;

@{ $data_ref->{'set'}->{'L3_indivs'} }  = (9..13);
@{ $data_ref->{'set'}->{'L4_indivs'} }  = (14..18);
@{ $data_ref->{'set'}->{'nhr_indivs'} } = (19..23);

my @sets = qw( L3_indivs L4_indivs nhr_indivs );

print "\n";

foreach my $column_set (@sets) { 
    foreach my $i (0..4) { 
        my $column_index = ${ $data_ref->{'set'}->{$column_set} }[$i];
        my $column_label = `cut -f $column_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | head --lines=1`;
        chomp $column_label;
        $column_label =~ s/_nz\z//;
        $column_label =~ s/indiv/i/;
        print "\t$column_label";
    }
    print "\n";

    foreach my $query_set (@sets) {
        foreach my $i (0..4) {
            my $query_index = ${ $data_ref->{'set'}->{$query_set} }[$i];
            my $query_label = `cut -f $query_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | head --lines=1`;
            chomp $query_label;
            my @query_data  = `cut -f $query_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | tail --lines=+2`;
            chomp @query_data;
            $query_label =~ s/_nz\z//; 
            $query_label =~ s/indiv/i/;
            print "$query_label";
            foreach my $j (0..4) { 
                my $column_index = ${ $data_ref->{'set'}->{$column_set} }[$j];
                my @column_data = `cut -f $column_index LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | tail --lines=+2`;
                chomp @column_data;

                my $lineFit = Statistics::LineFit->new();
                $lineFit->setData (\@query_data, \@column_data) or die "Invalid data";
                my ($intercept, $slope) = $lineFit->coefficients();
                defined $intercept or die "Can't fit line if x values are all equal";
                my $rSquared = $lineFit->rSquared();
                $rSquared = sprintf("%.2f", $rSquared);
                print "\t$rSquared";
            }
            print "\n";
        }
        print "\n";
    }
}
