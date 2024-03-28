#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;

my $data_ref;

my $header = "Seq\tMean\tMedian\tMax\tMin";

while ( my $input = <> ) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (\S+) \t \d+ \t \d+ \t (\S+) \z/xms ) ) {
        my $seq = $1;
        my $val = $2;
        if (! looks_like_number($val) ) {
            die "Value \"$val\" is not a number in: $input\n";
        }
        push @{ $data_ref->{'seq'}->{$seq}->{'val'} }, $val;
    }
    elsif ( $input !~ /\A [#] /xms ) {
       die "Cannot parse: $input\n";
    }
}

my @seqs = sort keys %{ $data_ref->{'seq'} };
foreach my $seq (@seqs) {
    my @vals = @{ $data_ref->{'seq'}->{$seq}->{'val'} };

    my $stat1 = Statistics::Descriptive::Full->new();
    $stat1->add_data(@vals);

    my $vals_min = $stat1->min();
    $vals_min    = commify($vals_min);

    my $vals_max = $stat1->max();
    $vals_max    = commify($vals_max);

    my $vals_mean = $stat1->mean();
    $vals_mean    = sprintf("%.1f", $vals_mean);
    $vals_mean    = commify($vals_mean);

    my $vals_median = $stat1->median();
    $vals_median    = sprintf("%.1f", $vals_median);
    $vals_median    = commify($vals_median);

    print "$header\n" if $header;
    $header = q{};

    print "$seq\t$vals_mean\t$vals_median\t$vals_max\t$vals_min\n";
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

