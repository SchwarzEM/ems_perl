#!/usr/bin/env perl

# count_ecoli_nt.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/4/2009.
# Purpose: get basic statistics from PS1010 assembly headers, in list or FASTA file.
# N.B.: totally depends on idiosyncracies of how the contigs are named.

use strict;
use warnings;
use Statistics::Descriptive;

my $nt_count;
my @sizes = ();

while (my $input = <>) { 
    if ( $input =~ / \A (?: > ){0,1} NODE_\d+_length_ (\d+) /xms ) { 
        $nt_count = $1;
        push @sizes, $nt_count;
    }
}

# Sort in ascending numerical order, so that homebrewed N50 subroutine can work.
@sizes = sort { $a <=> $b } @sizes;

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@sizes);

my $sum     = $stat->sum();    # Total nt of sequence.
my $rawsum = $sum;
$sum = commify($sum);

my $contigs = $stat->count();  # Total number of contigs.
$contigs = commify($contigs);

my $mean    = $stat->mean();
$mean       = sprintf("%.2f", $mean);

my $std_dev = $stat->standard_deviation();
$std_dev    = sprintf("%.2f", $std_dev);

my $min     = $stat->min();
$min        = commify($min);

my $max     = $stat->max();
$max        = commify($max);

my $median = $stat->median();
$median    = sprintf("%.2f", $median);

my $n50 = get_n50(\@sizes,$rawsum);
$n50    = sprintf("%.2f", $n50);

print "Total:  $sum nt in $contigs contigs.\n";
print "Mean:   $mean; std. dev. $std_dev; min. $min; max. $max\n";
print "Median: $median\n";
print "N50:    $n50\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

# Assumes array sorted in ascending values, with precomputed total value:
sub get_n50 { 
    my ($sizes_ref, $sum) = @_;
    my $midsum = ($sum / 2);
    my $tally = 0;
    my @_sizes = @{ $sizes_ref };
    foreach my $value (@_sizes) { 
        $tally += $value;
        if ($tally > $midsum) { 
            return $value;
        }
    }
    return;
}

