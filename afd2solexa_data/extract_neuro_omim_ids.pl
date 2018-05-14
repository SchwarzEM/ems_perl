#!/usr/bin/env perl

# extract_neuro_omim_ids.pl -- Erich Schwarz <emsch@caltech.edu>, 4/28/2012.
# Purpose: given OMIM neuro tables screen-scraped from omim.org, list their members as OMIM:\d+.

use strict;
use warnings;

my $omim_id  = q{};
my $omim_txt = q{};
my $extra    = q{};
my $reading  = 0;

while (my $input = <>) { 
    # 1 :     # 181500. SCHIZOPHRENIA; SCZD <
    if ( (! $reading) and ( $input =~ /\A \d /xms ) and ( $input =~ /\A \S*\d \s+ : \s .+? (\d+)\. \s (.+) \z /xms ) ) { 
        $omim_id  = $1;
        $omim_txt = $2;
        $omim_txt =~ s/\A\s+//;
        $omim_txt =~ s/\s+\z//;
        $reading  = 1;
        if ( $omim_txt =~ / \A (.+) \S{2} clinicalSynopsis /xms ) {
            $omim_txt = $1;
            $reading  = 0;
            $omim_txt =~ s/\s+\z//;
            print 'OMIM:', $omim_id, "\t$omim_txt\n";
            $omim_id  = q{};
            $omim_txt = q{};
            $extra    = q{};
        }
    }
    elsif ($reading) { 
        if ( $input =~ /\A (.*) \S{2} clinicalSynopsis /xms ) { 
            $extra    = $1;
            $reading  = 0;
            $omim_txt = "$omim_txt $extra";
            $omim_txt =~ s/\s+\z//;
            print 'OMIM:', $omim_id, "\t$omim_txt\n";
            $omim_id  = q{};
            $omim_txt = q{};
            $extra    = q{};
        }
        elsif ( $input =~ / \S /xms ) { 
            $input =~ s/\A\s+//;
            $input =~ s/\s+\z//; 
            $omim_txt = "$omim_txt $extra";
        }
        else { 
            die "Can't parse OMIM:$omim_id -- $omim_txt [and] $input\n";
        }
    } 
}

