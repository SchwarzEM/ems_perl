#!/usr/bin/env perl

use strict;
use autodie;
use warnings;

my $v1_to_v2 = q{};
my $v1_amps  = q{};

$v1_to_v2 = $ARGV[0] if $ARGV[0];
$v1_amps  = $ARGV[1] if $ARGV[1];

my $header = "Gene\tAMP";

my $data_ref;

if ( (! $v1_to_v2 ) or (! $v1_amps ) ) {
    die "Format: acey_amps_v1.to.v2.1_18nov2023.pl [v1 to v2.1 gene map] [v1 amp table] > [v2.1 genes => amps]\n";
}

open my $V1_TO_V2, '<', $v1_to_v2;
while ( my $input = <$V1_TO_V2> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S.*) \z/xms ) {
        my $v1      = $1;
        my $v2_list = $2;

        if ( exists $data_ref->{'v1'}->{$v1} ) {
            die "In v1_to_v2 file $v1_to_v2, redundant annotation of: $v1\n";
        }

        if ( $v1 ne 'Gene' ) {
            $data_ref->{'v1'}->{$v1}->{'v2_list'} = $v2_list;
        }
    }
    else {
        die "From v1_to_v2 file $v1_to_v2, cannot parse: $input\n";
    }
}
close $V1_TO_V2;

open my $V1_AMPS, '<', $v1_amps;
while ( my $input = <$V1_AMPS> ) {
    chomp $input;

    my $amp = q{};
    my $v1  = q{};

    if ( ( $input !~ / \A AMP [ ] Group /xms ) and ( $input =~ /\A (.+) \s+ \S+ \s+ (\S+) \z/xms ) ) {
        $amp = $1;
        $v1  = $2;

        if (! exists $data_ref->{'v1'}->{$v1}->{'v2_list'} ) {
            die "Cannot map v1 gene $v1 to v2.1 gene(s)\n";
        }

        my $v2_list = $data_ref->{'v1'}->{$v1}->{'v2_list'};
        my @v2s = split '; ', $v2_list;
        foreach my $v2 (@v2s) {
            print "$header\n" if $header;
            $header = q{};

            # Strip off any flanking whitespace that I may happen to be getting for the AMP annotation.
            $amp =~ s/\A\s+//;
            $amp =~ s/\s+\z//;

            print "$v2\t$amp\n";
        }
    }
    elsif ( $input !~ / \A AMP [ ] Group /xms ) {
        die "From v1_amps table $v1_amps, cannot parse: $input\n";
    }
}
close $V1_AMPS;

