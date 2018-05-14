#!/usr/bin/env perl

# summarize_nucmer_coords.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/7/2012; but very important debugging on 5/31/2016.
# Purpose: given a nucmer coords file, summarize -- for each query scaffold -- what its percentages of match are.  Mainly useful in assessing (in-)consistent chromosomal affinity.

use strict;
use warnings;

my $data_ref;

# Ignore the first five lines if they are *not* data.
my $lines = 0;

my $header = "Query\tSkew\tMatches";

while (my $input = <>) { 

# Example of first six lines of input:

# /mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/mummer/c_briggsae.PRJNA13758.WS254.genomic_masked.chr_only.sorted.fa /mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/mummer/nigoni_2015.12.01_gen.dna.fa
# NUCMER
# 
#     [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  | [TAGS]
# ==========================================================================================================
#     2007     2693  |  1678427  1677729  |      687      699  |    83.00  | 15455979  2935987  | Cbri_I	nigoni_2015.12.01_016

# Example of last three lines of input:

# 21531223 21534093  |   186153   183240  |     2871     2914  |    88.89  | 21540570  1366633  | Cbri_X	nigoni_2015.12.01_031
# 21531223 21534093  |    20012    22913  |     2871     2902  |    89.32  | 21540570    39279  | Cbri_X	nigoni_2015.12.01_097
# 21533428 21534080  |  2210922  2211574  |      653      653  |    89.28  | 21540570  2986942  | Cbri_X	nigoni_2015.12.01_015

    chomp $input;
    $lines++;
    if ( $lines >= 6 ) { 
        # Note that we start the line pattern after "\A" with "\s*", *not* "\s+"; failing to do this caused erroneous data loss in 2012.

        if ( $input =~ / \A \s* \d+ \s+ \d+   \s+  \|  # [S1]     [E1]
                            \s+ \d+ \s+ \d+   \s+  \|  # [S2]     [E2]
                            \s+ \d+ \s+ (\d+) \s+  \|  # [LEN 1]  [LEN 2]=$1
                            \s+ \d+\.\d+   \s+  \|  # [% IDY]
                         \s+ \d+ \s+ (\d+) \s+ \|   # [LEN R]  [LEN Q]=$2 
                         \s+ (\S+) \s+ (\S+) \s*    # [TAGS] = $3, $4
                         \z /xms ) {
            my ($len_match, $len_query, $match, $query) = ($1, $2, $3, $4);

            # Enforce unique, well-defined length for query sequence:
            if ( ( exists $data_ref->{'query'}->{$query}->{'len_query'} ) and ( $len_query != $data_ref->{'query'}->{$query}->{'len_query'} ) ) { 
                die "Can't parse: $input\n";
            }
            if (! exists $data_ref->{'query'}->{$query}->{'len_query'} ) { 
                $data_ref->{'query'}->{$query}->{'len_query'} = $len_query;
            } 

            # Length of matches to query will always add, in match-specific way:
            $data_ref->{'query'}->{$query}->{'match'}->{$match}->{'len_match'} += $len_match;
        }
        # Absolutely require parsing, thus ensuring that no data will be lost.
        else {
            die "Can't parse input line: $input\n";
        }
    }
}

foreach my $query ( sort keys %{ $data_ref->{'query'} } ) { 
    my $match_text = q{};
    my @match_summaries = ();

    my $len_query = $data_ref->{'query'}->{$query}->{'len_query'};
    # Archive this just once per query, so I can eventually summarize everything:
    $data_ref->{'total_query_lens'} += $len_query;

    my $match1to2_ratio = 'n/a';
    my $match1_nt       = 0;
    my $match2_nt       = 0;
    my @match_nt_vals   = ();

    if ( ( exists $data_ref->{'query'}->{$query}->{'match'} ) and ( exists $data_ref->{'query'}->{$query}->{'len_query'} ) ) {  
        my $i = 0;
        my @matches = sort {     $data_ref->{'query'}->{$query}->{'match'}->{$b}->{'len_match'} 
                             <=> $data_ref->{'query'}->{$query}->{'match'}->{$a}->{'len_match'} } 
                           keys %{ $data_ref->{'query'}->{$query}->{'match'} };
        foreach my $match (@matches) { 
            $i++;
            my $len_match = $data_ref->{'query'}->{$query}->{'match'}->{$match}->{'len_match'};

            # Archive this by position of match, so that we can eventually summarize the *total* bias of matches from best to worst, for the whole assembly:
            $data_ref->{'total_match_position'}->{$i}->{'total_match_len'} += $len_match;

            # Archive these by descending value
            push @match_nt_vals, $len_match;

            # Now get human-readable individual summaries:
            my $percent_match = ( $len_match * 100 / $len_query );
            my $len_match_p = commify($len_match);
            $percent_match = sprintf ( "%.2f", $percent_match );
            my $match_summary = "$match $percent_match" . '% [' . "$len_match_p nt]";
            push @match_summaries, $match_summary;
        }
        $match_text = join '; ', @match_summaries;
    }

    if ( exists $match_nt_vals[0] ) {
        $match1_nt = $match_nt_vals[0];
    }
    if ( exists $match_nt_vals[1] ) {
        $match2_nt = $match_nt_vals[1];
    }

    # If there is a real ratio, compute it; if it is effectively infinite, note that
    if ( ( $match1_nt >= 1 ) and ( $match2_nt == 0 ) ) {
        $match1to2_ratio = 'Inf';
    }
    if ( ( $match1_nt >= 1 ) and ( $match2_nt >= 1 ) ) {
        $match1to2_ratio = ($match1_nt/$match2_nt);
        $match1to2_ratio = sprintf ( "%.2f", $match1to2_ratio );
    }

    my $len_query2 = $data_ref->{'query'}->{$query}->{'len_query'};
    $len_query2 = commify($len_query2);

    print "$header\n" if $header;
    $header = q{};
    print "$query [$len_query2 nt]\t$match1to2_ratio\t$match_text\n";
}

print "\n";
my $total_nt = $data_ref->{'total_query_lens'};
my $total_nt_p = commify($total_nt);
print "Total seq.: $total_nt_p nt\n";

my @match_positions = sort { $a <=> $b } keys %{ $data_ref->{'total_match_position'} };
foreach my $pos (@match_positions) { 
    my $total_match_len = $data_ref->{'total_match_position'}->{$pos}->{'total_match_len'};
    my $total_match_perc = ( $total_match_len * 100 / $total_nt );
    $total_match_perc = sprintf ( "%.2f", $total_match_perc );
    my $total_match_len_p = commify($total_match_len);
    print "Match $pos: $total_match_perc% [$total_match_len_p nt]\n";
}
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

