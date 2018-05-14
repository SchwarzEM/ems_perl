#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::LineFit;

my $data_ref;

my @sets1 = (6..23);
my @sets2 = (2..5);
my @all_sets = @sets1;
push @all_sets, @sets2;

foreach my $set (@all_sets) { 
    $data_ref->{'set'}->{$set}->{'label'}     = `cut -f $set LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | head --lines=1`;
    @{ $data_ref->{'set'}->{$set}->{'data'} } = `cut -f $set LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | tail --lines=+2`;
    chomp $data_ref->{'set'}->{$set}->{'label'};
    chomp @{ $data_ref->{'set'}->{$set}->{'data'} };
}

print "\n";
foreach my $set1 (@sets1) { 
    foreach my $set2 (@sets2) {
        my $label1 = $data_ref->{'set'}->{$set1}->{'label'};
        my $label2 = $data_ref->{'set'}->{$set2}->{'label'};
        my @data1  = @{ $data_ref->{'set'}->{$set1}->{'data'} };
        my @data2  = @{ $data_ref->{'set'}->{$set2}->{'data'} };

        my $lineFit = Statistics::LineFit->new();
        $lineFit->setData (\@data1, \@data2) or die "Invalid data";
        my ($intercept, $slope) = $lineFit->coefficients();
        defined $intercept or die "Can't fit line if x values are all equal";
        my $rSquared = $lineFit->rSquared();
        $rSquared = sprintf("%.2f", $rSquared);

        print "$label1 vs. $label2, r^2:\t$rSquared\n";
    }
    print "\n";
}



