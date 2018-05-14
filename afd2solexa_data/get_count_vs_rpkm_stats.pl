#!/usr/bin/env perl

# get_count_vs_rpkm_stats.pl -- Erich Schwarz <emsch@caltech.edu>, 10/15/2011.
# Purpose: given a simple table of wordcount vs. rpkm values, make simple binned distributions of rpkm values for each value of wordcount.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;

my $data_ref; 

while (my $input = <>) { 
    chomp $input;
    my $wordcount = q{};
    my $rpkm      = q{};
    if ( $input =~ /\A (\d+) \t (\S+) \z /xms ) { 
        $wordcount = $1;
        $rpkm      = $2;
        if ( (! looks_like_number($wordcount) )
             or (! looks_like_number($rpkm)   ) ) { 
            die "Can't see numerical values in: $input\n";
        }
        push @{ $data_ref->{'wordcount'}->{$wordcount}->{'rpkm'} }, $rpkm;
    }
}

my @wordcounts1 = sort { $a <=> $b } keys %{ $data_ref->{'wordcount'} };

print "Words\tN\tMean\ts.d.\tMin.\tMedian\tMax.\n";

foreach my $wordcount1 (@wordcounts1) {  
    my @rpkms1 = @{ $data_ref->{'wordcount'}->{$wordcount1}->{'rpkm'} };

    my $stat1 = Statistics::Descriptive::Full->new();
    $stat1->add_data(@rpkms1);

    my $bin_n = $stat1->count();  # Total number of RPKM values observed for a given word count.
    $bin_n    = commify($bin_n);

    my $bin_mean    = $stat1->mean();
    $bin_mean       = sprintf("%.2f", $bin_mean);
    $bin_mean       = commify($bin_mean);

    my $bin_std_dev = $stat1->standard_deviation();
    $bin_std_dev    = sprintf("%.2f", $bin_std_dev);
    $bin_std_dev    = commify($bin_std_dev);

    my $bin_min     = $stat1->min();
    $bin_min        = commify($bin_min);

    my $bin_median = $stat1->median();
    $bin_median    = sprintf("%.2f", $bin_median);
    $bin_median    = commify($bin_median);

    my $bin_max     = $stat1->max();
    $bin_max        = commify($bin_max);

    print "$wordcount1\t$bin_n\t$bin_mean\t$bin_std_dev\t$bin_min\t$bin_median\t$bin_max\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


