#!/usr/bin/perl

# list_WB_macroannots.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/9/2007.
# Purpose: nonredudantly list macro-annot refs., and find which ones are WormBase papers.

# One example of usage:
#   ./list_WB_macroannots.pl   \
#       /home/schwarz/taygeta/science2/wormbase/Ongoing/wpa_xref.txt    \
#       /home/schwarz/taygeta/science2/macro_sets/working.macro-annot-list.txt > test.out

use strict;
use warnings;

my %pmid2wbpaper = ();
my %seen_pmid    = ();

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /(WBPaper\d+).+?(pmid\d+)/) { 
        $pmid2wbpaper{$2} = $1;
    }
    elsif ($input =~ /PMID:\s+(\d+)/) { 
       $seen_pmid{"pmid$1"} = 1;
    }
    elsif ($input =~ /(pmid\d+)/) {
       $seen_pmid{"$1"} = 1;
    }
}

foreach my $pmid (sort keys %seen_pmid) { 
    if ($pmid2wbpaper{$pmid}) { 
        print "$pmid\t$pmid2wbpaper{$pmid}\n";
    }
    else { 
        print "$pmid\n"
    }
}

