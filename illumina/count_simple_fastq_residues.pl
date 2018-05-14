#!/usr/bin/env perl

# count_simple_fastq_residues.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/28/2010. 
# Purpose: get basic statistics from simple FASTQ files; counts residues directly.

use strict;
use warnings;
use Statistics::Descriptive;

my $nt_count = 0;
my @sizes = ();

my $i = 0;

while (my $input = <>) { 
    chomp $input;
    $i++;
    # Use !(modulus) function to get sequence from every 2cd line of the 4-line records.
    if ( ( ( $i % 4 ) == 2 ) and ( $input =~ / \A \s* [ACGTNacgt]+ \s* \z /xms ) ) { 
        $input =~ s/\s//g;
        $nt_count += length($input);
        push @sizes, $nt_count;
        $nt_count = 0;
    }
    if ( ( ( $i % 4 ) == 2 ) and ( $input !~ / \A \s* [ACGTNacgt]+ \s* \z /xms ) ) {
        die "Unexpected sequence: $input\n";
    }
}

# Sort in ascending numerical order, so that homebrewed N50 subroutine can work.
@sizes = sort { $a <=> $b } @sizes;

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@sizes);

my $sum     = $stat->sum();    # Total nt of sequence.
my $rawsum = $sum;
$sum = commify($sum);

my $reads = $stat->count();  # Total number of reads.
$reads = commify($reads);

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

print "Total:  $sum nt in $reads reads.\n";
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

