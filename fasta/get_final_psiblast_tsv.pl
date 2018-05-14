#!/usr/bin/env perl

# get_final_psiblast_tsv.pl -- Erich Schwarz <ems394@cornell.edu>, 5/3/2014.
# Purpose: from stock tabular '-m 6' psiblast output, get the *last* round of hits; revised to deal with domain/motif searches on proteomes (which can hit query sequence 2+ times).

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my $query_seq    = q{};
my $target_seq   = q{};
my $percent_hit  = 0;
my @result_lines = ();

while (my $input = <>) { 
    # Sample input line for query sequence re-hitting itself: 
    # Acey_s0004.g1962.t1     Acey_s0004.g1962.t1     100.00	[etc.]
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \t .+ \z/xms ) { 
        $query_seq          = $1;
        $target_seq         = $2;
        my $new_percent_hit = $3;

        # Enforce correct number:
        if ( (! looks_like_number($new_percent_hit) ) or ( $new_percent_hit < 0 ) or ( $new_percent_hit > 100 ) ) { 
            die "Can't parse % hit value ($new_percent_hit) in input line: $input\n";
        }

        # This should empirically force the $percent_hit to always be the *largest* number possible.
        # If it is somehow set to anything less than 100, the script issues a warning (*every* time that happens).
        # The standard, non-pathological thing that should happen is that $percent_hit should become 100.00 on the very first input line,
        #     and then stay there from then on.

        if ( $new_percent_hit > $percent_hit ) {
            $percent_hit = $new_percent_hit;
            if ( $percent_hit != 100 ) { 
                warn "Set %-hit threshold for start of new result stanza to $percent_hit (i.e. something other than 100) in: $input\n";
            }
        }
        if ( ( $query_seq eq $target_seq ) and ( $percent_hit = $new_percent_hit ) )  { 
            @result_lines = ();
        }
        push @result_lines, $input;
    }
    elsif ( $input =~ /\ASearch [ ] has [ ] CONVERGED!/xms ) { 
        warn $input;
    }
}

print @result_lines;

