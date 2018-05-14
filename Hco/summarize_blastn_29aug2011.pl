#!/usr/bin/env perl

# summarize_blastn_29aug2011.pl -- Erich Schwarz <emsch@caltech.edu>, 8/29/2011.
# Purpose: given a set of similarities via BlastN, compute the overall similarity for one sequence to another.

use strict;
use warnings;

my $seqfile   = q{};
my $blastfile = q{};

$seqfile   = $ARGV[0] if $ARGV[0];
$blastfile = $ARGV[1] if $ARGV[1];

my $seq1       = q{};
my $seq2       = q{};
my $diffs      = q{};
my $seq1_start = q{};
my $seq1_end   = q{};
my $seq2_start = q{};
my $seq2_end   = q{};

my $data_ref;

open my $SEQ, '<', $seqfile or die "Can't open sequence file $seqfile: $!";
while (my $input = <$SEQ>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $seq1 = $1;
    }
    elsif ( $input =~ / [ACGTNactgn ] /xms ) { 
        $input =~ s/[^ACGTNactgn]//g;
        my $seq_len = length($input);
        $data_ref->{'seq'}->{$seq1}->{'length'} += $seq_len;
    }
}
close $SEQ or die "Can't close handle to sequence file $seqfile: $!";

# $1              $2                              $3              $4      $5      $6      $7
# Hco_McM_000001  Hco_McM_002363  81.41   1958    241     17      769756  771665  59626   61508   0.0     1925
# Hco_McM_000001  Hco_McM_002363  79.28   2013    245     29      779206  781105  74679   76632   0.0     1775

open my $BLAST, '<', $blastfile or die "Can't open BlastN file $blastfile: $!";
while (my $input = <$BLAST>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t 
                       \d{2,3}\.\d{2} \t \d+ \t 
                       (\d+) \t 
                       \d+ \t 
                       (\d+) \t (\d+) \t 
                       (\d+) \t (\d+) \t /xms ) { 
        $seq1       = $1;
        $seq2       = $2;
        $diffs      = $3;
        $seq1_start = $4;
        $seq1_end   = $5;
        $seq2_start = $6;
        $seq2_end   = $7;
        my $seq_pair = "$seq1\t$seq2";
        $data_ref->{'seq_pair'}->{$seq_pair}->{'diffs'} += $diffs;

        if ( $seq1_end < $seq1_start ) { 
            ($seq1_start, $seq1_end) = ($seq1_end, $seq1_start);
        }
        if ( $seq2_end < $seq2_start ) {
            ($seq2_start, $seq2_end) = ($seq2_end, $seq2_start);  
        }

        foreach my $i ($seq1_start..$seq1_end) { 
            $data_ref->{'seq_pair'}->{$seq_pair}->{'seq1'}->{$seq1}->{'match'}->{$i} = 1;
        }
        foreach my $i ($seq2_start..$seq2_end) {
            $data_ref->{'seq_pair'}->{$seq_pair}->{'seq2'}->{$seq2}->{'match'}->{$i} = 1;
        }
    }
    else { 
        die "Can't parse: $input\n";
    }
}
close $BLAST or die "Can't close handle to BlastN file $blastfile: $!";

my $header = "Seq1\tSeq1_matches\tSeq1_len\tSeq1_%match\tSeq2\tSeq2_matches\tSeq2_len\tSeq2_%match\tDiffs\tSeq2_%diffs/matches";
print "$header\n";

foreach my $seq_pair1 (sort keys %{ $data_ref->{'seq_pair'} } ) { 
    if ( $seq_pair1 !~ /\A \S+ \t \S+ \z/xms ) { 
        die "Can't parse sequence pair: $seq_pair1\n";
    }
    if ( $seq_pair1 =~ /\A (\S+) \t (\S+) \z/xms ) {
        $seq1 = $1;
        $seq2 = $2;
    }
    my $seq1_length = $data_ref->{'seq'}->{$seq1}->{'length'};
    my $seq2_length = $data_ref->{'seq'}->{$seq2}->{'length'};

    my $seq1_matchcount = keys %{ $data_ref->{'seq_pair'}->{$seq_pair1}->{'seq1'}->{$seq1}->{'match'} };
    my $seq2_matchcount = keys %{ $data_ref->{'seq_pair'}->{$seq_pair1}->{'seq2'}->{$seq2}->{'match'} };

    my $seq1_frac_match = (100 * $seq1_matchcount / $seq1_length);
    $seq1_frac_match = sprintf("%.2f", $seq1_frac_match);    

    my $seq2_frac_match = (100 * $seq2_matchcount / $seq2_length);
    $seq2_frac_match = sprintf("%.2f", $seq2_frac_match);

    $diffs = 0;
    if ( exists $data_ref->{'seq_pair'}->{$seq_pair1}->{'diffs'} ) { 
        $diffs = $data_ref->{'seq_pair'}->{$seq_pair1}->{'diffs'};
    }
    my $seq2_frac_diffs = (100 * $diffs / $seq2_matchcount);
    $seq2_frac_diffs = sprintf("%.2f", $seq2_frac_diffs);

    $seq1_length     = commify($seq1_length);
    $seq1_matchcount = commify($seq1_matchcount);
    $seq2_length     = commify($seq2_length);
    $seq2_matchcount = commify($seq2_matchcount);
    $diffs           = commify($diffs);

    print "$seq1\t$seq1_matchcount\t$seq1_length\t$seq1_frac_match\t$seq2\t$seq2_matchcount\t$seq2_length\t$seq2_frac_match\t$diffs\t$seq2_frac_diffs\n";
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

