#!/usr/bin/env perl

# bowtie2size_dist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/10/2011.
# Purpose: given a bowtie-produced SAM alignment of paired-end reads, compute the mean, s.d., min., and max. of their distances.

use strict;
use warnings;
use Statistics::Descriptive;

my $pseq = q{};
my $dist = q{};
my $read = q{};

my @dists = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \S+ \t (?: 99|163) \t \S+ \t \d+ \t \d+ \t \S+ \t (\S+) \t \d+ \t (\d+) \t ([ACGTNacgtn]+) \t /xms ) { 
        $pseq = $1;
        $dist = $2;
        $read = $3;
        if ( ( $pseq ne '=' ) or ( $dist < 1 ) or ( $dist != int($dist) ) ) {  
            die "Can't parse input line: $input\n";
        }
        my $read_len = length($read);
        $dist += $read_len;
        push @dists, $dist;
    }
}

my $stat1 = Statistics::Descriptive::Full->new();
$stat1->add_data(@dists);

my $inserts = $stat1->count();  
$inserts    = commify($inserts);

my $dist_mean    = $stat1->mean();
$dist_mean       = sprintf("%.1f", $dist_mean);
$dist_mean       = commify($dist_mean);

my $dist_std_dev = $stat1->standard_deviation();
$dist_std_dev    = sprintf("%.1f", $dist_std_dev);
$dist_std_dev    = commify($dist_std_dev);

my $dist_min     = $stat1->min();
$dist_min        = commify($dist_min);

my $dist_max     = $stat1->max();
$dist_max        = commify($dist_max);

print "\n";
print "Inserts:        $inserts\n";
print "Mean len. nt:   $dist_mean\n";
print "Std. dev. nt:   $dist_std_dev\n";
print "Min. len. nt:   $dist_min\n";
print "Max. len. nt:   $dist_max\n";
print "\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

