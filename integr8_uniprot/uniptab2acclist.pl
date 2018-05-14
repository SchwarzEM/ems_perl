#!/usr/bin/env perl

# uniptab2acclist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/12/2008.
# Purpose: given 'UniProtKB AC/ID' to 'RefSeq' IDs from beta.uniprot.org, get sensible list for Batch Entrez.

use strict;
use warnings;

my %up2acc = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ / \A (\w+_\w+) \s+ (\w+_\w+?\S*) /xms ) { 
        my ($unipr_id, $acc) = ($1, $2); 

        # Anything is better than nothing:
        if (! exists $up2acc{$unipr_id} ) { 
            $up2acc{$unipr_id} = $acc;
        }

        # But it's better to settle on a real RefSeq ID (de-suffixed):
        if ($acc =~ / \w+\.\w+ /xms) { 
            $acc =~ s/\.\w+//;
            $up2acc{$unipr_id} = $acc;
        }
    }
}

foreach my $acc (sort values %up2acc) { 
    print "$acc\n";
}

