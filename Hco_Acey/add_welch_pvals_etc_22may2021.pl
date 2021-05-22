#!/usr/bin/env perl

# Erich Schwarz <ems394@cornell.edu>, 5/22/2021.

# Subroutines in this Perl script for computing Welch test p-values were copied from:
# "Using Burkhardt's 'incomplete beta'" in https://rosettacode.org/wiki/Welch%27s_t-test#Perl 

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;
use List::Util 'sum';

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

my $header = "Diff_all"
             . "\t"
             . "p-value_all"
             . "\t"
             . "Mean_all"
             . "\t"
             . "Max_all"
             . "\t"
             . "Min_all"
             . "\t"
             . "Diff_hw"
             . "\t"
             . "p-value_hw"
             . "\t"
             . "Mean_hw"
             . "\t"
             . "Max_hw"
             . "\t"
             . "Min_hw"
             . "\t"
             . "Mean_control"
             . "\t"
             . "Max_control"
             . "\t" 
             . "Min_control";

if (! $infile ) {
    die "Format: add_welch_pvals_19may2021.pl [infile] > [outfile, with added stats for each line]\n";
}

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A ( \S+ \t [^\t]+ \t \S+) \t  # Allow 2+ homology annots in column 2
                       ( (?:\S+\t){8}  \S+) \t
                       ( (?:\S+\t){23} \S+) \t
                       ( (?:\S+\t){4}  \S+)
                    \z/xms ) { 

        my $probe_data = $1;

        my $set_a_text = $2;   # all hookworm-infected data
        my $set_c_text = $3;   # control data
        my $precalcs   = $4;

        my @set_a = split /\t/, $set_a_text;
        my @set_c = split /\t/, $set_c_text;

        # Out of 9 hookworm-infection data, extract the last 3 (which are hookworm-only):
        my @set_h = ($set_a[6], $set_a[7], $set_a[8]);

        # Emit exactly one expanded header line at the beginning.
        if ( $header ) {
            $header = "$probe_data\t$set_a_text\t$set_c_text\t$header\t$precalcs";
            print "$header\n";
            $header = q{};
        }

        # After that, print out expanded data lines if possible.
        elsif ( looks_like_number($set_a[0]) ) {
            my $stat_a = Statistics::Descriptive::Full->new();
            my $stat_h = Statistics::Descriptive::Full->new();
            my $stat_c = Statistics::Descriptive::Full->new();

            $stat_a->add_data(@set_a);
            $stat_h->add_data(@set_h);
            $stat_c->add_data(@set_c);

            my $mean_a = $stat_a->mean();
            my $mean_h = $stat_h->mean();
            my $mean_c = $stat_c->mean();

            my $max_a = $stat_a->max();
            my $max_h = $stat_h->max();
            my $max_c = $stat_c->max();

            my $min_a = $stat_a->min();
            my $min_h = $stat_h->min();
            my $min_c = $stat_c->min();

            my $diff_ac = $mean_a - $mean_c ;
            my $diff_hc = $mean_h - $mean_c ;

            my $pvalue_a = calculate_P_value(\@set_a, \@set_c);
            my $pvalue_h = calculate_P_value(\@set_h, \@set_c);

            my $new_stats = "$diff_ac"
                            . "\t"
                            . "$pvalue_a"
                            . "\t"

                            . "$mean_a"
                            . "\t"
                            . "$max_a"
                            . "\t"
                            . "$min_a"
                            . "\t"

                            . "$diff_hc"
                            . "\t"
                            . "$pvalue_h"
                            . "\t"

                            . "$mean_h"
                            . "\t"
                            . "$max_h"
                            . "\t"
                            . "$min_h"
                            . "\t"

                            . "$mean_c"
                            . "\t"
                            . "$max_c"
                            . "\t"
                            . "$min_c"
                            ;

            print "$probe_data\t$set_a_text\t$set_c_text\t$new_stats\t$precalcs\n";
        }

        else {
            die "In input file $infile, can parse but cannot analyze: $input\n";
        }
    }
    else {
        die "In input file $infile, cannot parse: $input\n";
    }
}
close $INFILE;

# The following subroutines for computing Welch test p-values were copied from:
# "Using Burkhardt's 'incomplete beta'" in https://rosettacode.org/wiki/Welch%27s_t-test#Perl 

sub lgamma {
  my $x = shift;
  my $log_sqrt_two_pi = 0.91893853320467274178;
  my @lanczos_coef = (
      0.99999999999980993, 676.5203681218851, -1259.1392167224028,
      771.32342877765313, -176.61502916214059, 12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 );
  my $base = $x + 7.5;
  my $sum = 0;
  $sum += $lanczos_coef[$_] / ($x + $_)  for reverse (1..8);
  $sum += $lanczos_coef[0];
  $sum = $log_sqrt_two_pi + log($sum/$x) + ( ($x+0.5)*log($base) - $base );
  $sum;
}

sub calculate_P_value {
    my ($array1,$array2) = (shift, shift);
    return 1 if @$array1 <= 1 or @$array2 <= 1;
 
    my $mean1 = sum(@$array1);
    my $mean2 = sum(@$array2);
    $mean1 /= scalar @$array1;
    $mean2 /= scalar @$array2;
    return 1 if $mean1 == $mean2;
 
    my ($variance1,$variance2);
    $variance1 += ($mean1-$_)**2 for @$array1;
    $variance2 += ($mean2-$_)**2 for @$array2;
    return 1 if $variance1 == 0 and $variance2 == 0;
 
    $variance1 = $variance1/(@$array1-1);
    $variance2 = $variance2/(@$array2-1);
    my $Welch_t_statistic = ($mean1-$mean2)/sqrt($variance1/@$array1+$variance2/@$array2);
    my $DoF = (($variance1/@$array1+$variance2/@$array2)**2) /
                (
                 ($variance1*$variance1)/(@$array1*@$array1*(@$array1-1)) +
                 ($variance2*$variance2)/(@$array2*@$array2*(@$array2-1))
                );
    my $A     = $DoF / 2;
    my $value = $DoF / ($Welch_t_statistic**2 + $DoF);
    return $value if $A <= 0 or $value <= 0 or 1 <= $value;
 
    # from here, translation of John Burkhardt's C code
    my $beta = lgamma($A) + 0.57236494292470009 - lgamma($A+0.5); # constant is lgamma(.5), but more precise than 'lgamma' routine
    my $eps = 10**-15;
    my($ai,$cx,$indx,$ns,$pp,$psq,$qq,$qq_ai,$rx,$term,$xx);
 
    $psq = $A + 0.5;
    $cx = 1 - $value;
    if ($A < $psq * $value) { ($xx, $cx, $pp, $qq, $indx) = ($cx,   $value, 0.5,  $A, 1) }
    else                    { ($xx,      $pp, $qq, $indx) = ($value,         $A, 0.5, 0) }
    $term = $ai = $value = 1;
    $ns = int $qq + $cx * $psq;
 
    # Soper reduction formula
    $qq_ai = $qq - $ai;
    $rx = $ns == 0 ? $xx : $xx / $cx;
    while (1) {
        $term = $term * $qq_ai * $rx / ( $pp + $ai );
        $value = $value + $term;
        $qq_ai = abs($term);
        if ($qq_ai <= $eps && $qq_ai <= $eps * $value) {
           $value = $value * exp ($pp * log($xx) + ($qq - 1) * log($cx) - $beta) / $pp;
           $value = 1 - $value if $indx;
           last;
        }
        $ai++;
        $ns--;
        if ($ns >= 0) {
            $qq_ai = $qq - $ai;
            $rx = $xx if $ns == 0;
        } else {
            $qq_ai = $psq;
            $psq = $psq + 1;
        }
    }
    $value
}

