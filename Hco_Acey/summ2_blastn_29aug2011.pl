#!/usr/bin/env perl

# summ2_blastn_29aug2011.pl -- Erich Schwarz <emsch@caltech.edu>, 8/29/2011.
# Purpose: given a set of similarities via BlastN, compute the fraction of each sequence that has matches to any other.

use strict;
use warnings;
use Getopt::Long;

my @seqfiles  = ();
my $seqfile   = q{};
my $blastfile = q{};
my $help;

GetOptions ( 'seqs=s{,}' => \@seqfiles,
             'blast=s'   => \$blastfile,
             'help'      => \$help, );

if ( $help or (! @seqfiles) or (! $blastfile) ) { 
    die "Format: summ2_blastn_29aug2011.pl --seqs|-s [1+ FASTA seqs.] --blast|-b [tabular BlastN] --help|-h [print this message]\n";
}


my $seq1       = q{};
my $seq2       = q{};
my $seq1_start = q{};
my $seq1_end   = q{};
my $seq2_start = q{};
my $seq2_end   = q{};

my $data_ref;

foreach my $seqfile (@seqfiles) { 
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
}

# $1              $2                                              $3      $4      $5      $6
# Hco_McM_000001  Hco_McM_002363  81.41   1958    241     17      769756  771665  59626   61508   0.0     1925
# Hco_McM_000001  Hco_McM_002363  79.28   2013    245     29      779206  781105  74679   76632   0.0     1775

open my $BLAST, '<', $blastfile or die "Can't open BlastN file $blastfile: $!";
while (my $input = <$BLAST>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t 
                       \d{2,3}\.\d{2} \t \d+ \t \d+ \t \d+ \t 
                       (\d+) \t (\d+) \t 
                       (\d+) \t (\d+) \t /xms ) { 
        $seq1       = $1;
        $seq2       = $2;
        $seq1_start = $3;
        $seq1_end   = $4;
        $seq2_start = $5;
        $seq2_end   = $6;

        if ( $seq1_end < $seq1_start ) { 
            ($seq1_start, $seq1_end) = ($seq1_end, $seq1_start);
        }
        if ( $seq2_end < $seq2_start ) {
            ($seq2_start, $seq2_end) = ($seq2_end, $seq2_start);  
        }

        foreach my $i ($seq1_start..$seq1_end) { 
            $data_ref->{'seq'}->{$seq1}->{'match'}->{$i} = 1;
        }
        foreach my $i ($seq2_start..$seq2_end) {
            $data_ref->{'seq'}->{$seq2}->{'match'}->{$i} = 1;
        }
    }
    else { 
        die "Can't parse: $input\n";
    }
}
close $BLAST or die "Can't close handle to BlastN file $blastfile: $!";

my $header = "Seq1\tSeq_matches\tSeq_len\tSeq_%match\tNonred_len\tTotal_len";
print "$header\n";

my $total_len = 0;
my $nr_len    = 0;

foreach my $sequence1 (sort keys %{ $data_ref->{'seq'} } ) { 
    my $seq1_length = $data_ref->{'seq'}->{$sequence1}->{'length'};
    my $seq1_matchcount = keys %{ $data_ref->{'seq'}->{$sequence1}->{'match'} };
    my $seq1_frac_match = (100 * $seq1_matchcount / $seq1_length);
    $seq1_frac_match = sprintf("%.2f", $seq1_frac_match);    

    my $non_redundant_len = $seq1_length - $seq1_matchcount;

    $total_len += $seq1_length;
    $nr_len    += $non_redundant_len;

    $seq1_length     = commify($seq1_length);
    $seq1_matchcount = commify($seq1_matchcount);

    my $readable_total_len = commify($total_len);
    my $readable_nr_len    = commify($nr_len);

    print "$sequence1\t$seq1_matchcount\t$seq1_length\t$seq1_frac_match\t$readable_nr_len\t$readable_total_len\n";
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

