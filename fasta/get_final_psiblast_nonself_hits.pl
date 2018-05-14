#!/usr/bin/env perl

# get_final_psiblast_nonself_hits.pl -- Erich Schwarz <ems394@cornell.edu>, 12/31/2016.
# Purpose: from one (or more!) streamed annotated tabular '-m 7' psiblast output, get the last round of non-self hits *if* those hits are converged.

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $query_seq    = q{};
my $target_seq   = q{};
my @non_self_hits = ();

while (my $input = <>) { 
    chomp $input;

    # Sample input lines:
    # [...]
    # # Iteration: 5
    # FL81_17794-RA   Csp5_scaffold_naive_04733.g27387.t1     [etc.]
    # FL81_17794-RA   FL81_17794-RA   [etc.]
    # [...]
    # Search has CONVERGED!

    # Start each round with an empty hit list:
    if ( $input =~ /\A [#] [ ] Iteration: [ ] \d+ \s* \z/xms ) {
        @non_self_hits = ();
    }

    # Add non-self hits:
    elsif ( $input =~ /\A (\S+) \t (\S+) \t .+ \z/xms ) { 
        $query_seq          = $1;
        $target_seq         = $2;
        if ( $target_seq ne $query_seq ) {
            push @non_self_hits, $target_seq;
        }
    }

    # If and when convergence occurs, print a sorted non-redundant non-self hit list, then clear it
    #     (allowing 2+ input files to be dealt with in one stream):
    elsif ( $input =~ /\A Search [ ] has [ ] CONVERGED! \s* \z/xms ) { 
        @non_self_hits = sort @non_self_hits;
        @non_self_hits = uniq(@non_self_hits);
        foreach my $non_self_hit (@non_self_hits) {
            print "$non_self_hit\n";
        }
        @non_self_hits = ();
    }
}

