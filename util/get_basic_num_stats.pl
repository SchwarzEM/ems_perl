#!/usr/bin/env perl

# get_basic_num_stats.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/20/2011.
# Purpose: given a stream of numbers, compute: N, mean, median, min, max, s.d.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;

my @input_numbers = ();

while (my $input = <>) { 
    chomp $input;
    if ( looks_like_number($input) ) { 
        push @input_numbers, $input;
    }
    else { 
        die "$input does not look like a number.\n";
    }
}

my $stat1 = Statistics::Descriptive::Full->new();
$stat1->add_data(@input_numbers);

my $input_count = $stat1->count();  # Total number of scaffolds.
$input_count      = commify($input_count);

my $input_mean    = $stat1->mean();
$input_mean       = commify($input_mean);

my $input_median = $stat1->median();
$input_median    = commify($input_median);

my $input_min     = $stat1->min();
$input_min        = commify($input_min);

my $input_max     = $stat1->max();
$input_max        = commify($input_max);

my $input_std_dev = $stat1->standard_deviation();
$input_std_dev    = commify($input_std_dev);

print "N:       $input_count\n",
      "Mean:    $input_mean\n",
      "Median:  $input_median\n",
      "Min.:    $input_min\n",
      "Max.:    $input_max\n",
      "S.d.:    $input_std_dev\n",
      ;

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

