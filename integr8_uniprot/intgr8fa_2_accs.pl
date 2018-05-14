#!/usr/bin/env perl

# intgr8fa_2_accs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/12/2008.
# Purpose: given an Integr8-derived FASTA file, extract list of acc. nos. usable by NCBI's Batch (protein) Entrez.

# N.B. Assumes, but does not totally require, this sort of header line:
# >Q8FYZ9_BRUSU   Q8FYZ9   Q8FYZ9_BRUSU   BR1706   Brucella suis [NCBI_TaxID=29461]   100.B_suis_1.dat
# produced from a *.dat file by uniprot_to_clean_tfa.pl.

use strict;
use warnings;

my %usable_IDs = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A > (\w+)_\w+ \s+ (\w+)? \s* /xms ) { 
        my ($intgr8_id_stem, $acc_id) = ($1, $2);
        if ( ($intgr8_id_stem ne $acc_id) and ($acc_id =~ /[A-Z]/xms) ) {
            warn "Discordant inferred acc. nos.: $intgr8_id_stem vs. $acc_id, from $input\n";
        }
        if ($acc_id =~ /[A-Z]/xms) { 
            $usable_IDs{$acc_id} = 1;
        }
        if ($acc_id !~ /[A-Z]/xms) {
            $usable_IDs{$intgr8_id_stem} = 1;
        }
    }
}

foreach my $usable_ID (sort keys %usable_IDs) { 
    print "$usable_ID\n";
}

