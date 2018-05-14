#!/usr/bin/env perl

# get_final_psiblast_tsv_26mar2013.pl -- Erich Schwarz <ems394@cornell.edu>, 3/26/2013 -- LEGACY version, kept for reproducibility of previous work.
# Purpose: from stock tabular '-m 6' psiblast output, get the *last* round of hits.  Note that this only works reliably when the query sequence is a 
#    full-length member of the database being searched (so that it will have exactly one hit at 100% identity, time after time).
#    It does *not* perform reliabily if a database of sequences is being searched with a subsequence (e.g., a repeated motif or domain) of
#    one full-length sequence in that database, if the full-length sequence has 2x copies of the motif/domain.
#    This legacy version has therefore been superceded by a revised get_final_psiblast_tsv.pl

use strict;
use warnings;

my $query_seq    = q{};
my $target_seq   = q{};
my @result_lines = ();

while (my $input = <>) { 
    if ( $input =~ /\A (\S+) \t (\S+) \t .+ \z/xms ) { 
        $query_seq  = $1;
        $target_seq = $2;
        if ( $query_seq eq $target_seq ) { 
            @result_lines = ();
        }
        push @result_lines, $input;
    }
    elsif ( $input =~ /\ASearch [ ] has [ ] CONVERGED!/xms ) { 
        warn $input;
    }
}

print @result_lines;

