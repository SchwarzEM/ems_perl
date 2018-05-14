#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::LineFit;

my $data_ref;

my @sets1 = (2, 6, 9..13, 3, 7, 14..18, 4, 8, 19..23);
my @sets2 = (2, 6, 3, 7, 4, 8, 5);

my @all_sets = @sets1;
push @all_sets, @sets2;

foreach my $set (@all_sets) { 
    $data_ref->{'set'}->{$set}->{'label'}     = `cut -f $set LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | head --lines=1`;
    @{ $data_ref->{'set'}->{$set}->{'data'} } = `cut -f $set LC_table_23jun2012_10K_nz.pools_means_indivs.tsv.txt | tail --lines=+2`;
    chomp $data_ref->{'set'}->{$set}->{'label'};
    chomp @{ $data_ref->{'set'}->{$set}->{'data'} };
}

print "\n";

foreach my $set2 (@sets2) { 
    my $label2 = $data_ref->{'set'}->{$set2}->{'label'};
    $label2 =~ s/_nz\z//;
    $label2 =~ s/_/ /g;
    $label2 =~ s/indiv Mean/iM/;
    $label2 =~ s/pool/pl/;
    print "\t$label2";
}
print "\n";

foreach my $set1 (@sets1) { 
    my $label1 = $data_ref->{'set'}->{$set1}->{'label'};
    $label1 =~ s/_nz\z//;
    $label1 =~ s/_/ /g;
    print "$label1";
    foreach my $set2 (@sets2) {
        my @data1  = @{ $data_ref->{'set'}->{$set1}->{'data'} };
        my @data2  = @{ $data_ref->{'set'}->{$set2}->{'data'} };

        my $lineFit = Statistics::LineFit->new();
        $lineFit->setData (\@data1, \@data2) or die "Invalid data";
        my ($intercept, $slope) = $lineFit->coefficients();
        defined $intercept or die "Can't fit line if x values are all equal";
        my $rSquared = $lineFit->rSquared();
        $rSquared = sprintf("%.2f", $rSquared);

        print "\t$rSquared";
    }
    print "\n";
}

print "\n";

