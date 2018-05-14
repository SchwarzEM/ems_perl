#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw{uniq};

my $data_ref;
my $index   = q{};
my @matches = ();

# Sample input:
#
# >Cluster 0
# 0       32757aa, >A5X6X5... *
# >Cluster 1
# 0       31319aa, >F1R7N8... *
# 1       31319aa, >ENSDARP00000099532... at 100.00%

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > Cluster \s+ \d+ \s* \z/xms ) {
        &compile_and_clear_seqs();
    }
    elsif ( $input =~ /\A \d+ \s+ \d+ aa, \s+ > (\S+) \.\.\. \s+ [*] \s* \z/xms ) {
        $index = $1;
    }
    elsif ( $input =~ /\A \d+ \s+ \d+ aa, \s+ > (\S+) \.\.\. \s+ at \s+ \d+\.\d+ [%] \s* \z/xms ) {
        my $match = $1;
        push @matches, $match;
    }
}
# One last clearing of data after reading the whole file:
&compile_and_clear_seqs();

my @index_seqs = sort keys %{ $data_ref->{'index'} };
foreach my $index_seq (@index_seqs) {
    my @match_seqs = @{ $data_ref->{'index'}->{$index_seq}->{'matches'} };
    foreach my $match_seq (@match_seqs) {
        print "$match_seq\t$index_seq\n";
    }
}

sub compile_and_clear_seqs {
    if ( $index and (@matches) ) {
        @matches = sort @matches;
        @matches = uniq @matches;
        @{ $data_ref->{'index'}->{$index}->{'matches'} } = @matches;
    }
    $index   = q{};
    @matches = ();
}

