#!/usr/bin/env perl

# eggnog3_to_wormhum_orths.pl -- Erich Schwarz <emsch@caltech.edu>, 4/28/2012.
# Purpose: extract simple worm-human ortholog tables from meNOG.members.txt (obtained from http://eggnog.embl.de/version_3.0/data/downloads/meNOG.members.txt.gz).

# Note species code, from http://eggnog.embl.de/version_3.0/data/downloads/species.v3.txt:
# 6239 == elegans
# 9606 == sapiens

use strict; 
use warnings;

my $ogrp    = q{};
my $species = q{};
my $prot    = q{};

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (meNOG\d+) \s+ (\d+)\.(\S+) \s /xms ) { 
        $ogrp    = $1;
        $species = $2;
        $prot    = $3;
        if ( $species == 6239 ) {
            # Remove superfluous transcript isoform suffixes, to get get CDS gene names:
            if ( $prot =~ /\A (\S+\.\d+) [a-z]{0,1}\.\d+ \z/xms ) {
                $prot = $1;
            }
            # Remove isoform a-z suffixes from proteins, to get CDS gene names: 
            $prot =~ s/[a-z]\z//;

            $data_ref->{'ogrp'}->{$ogrp}->{'worm'}->{$prot} = 1;
        }
        elsif ( $species == 9606 ) { 
            $data_ref->{'ogrp'}->{$ogrp}->{'human'}->{$prot} = 1;
        }
    }
}

foreach my $ogrp1 ( sort keys %{ $data_ref->{'ogrp'} } ) { 
    if ( ( exists $data_ref->{'ogrp'}->{$ogrp1}->{'worm'} ) and ( exists $data_ref->{'ogrp'}->{$ogrp1}->{'human'} ) ) { 
        foreach my $wpep ( sort keys %{ $data_ref->{'ogrp'}->{$ogrp1}->{'worm'} } ) { 
            foreach my $hpep ( sort keys %{ $data_ref->{'ogrp'}->{$ogrp1}->{'human'} } ) {
                print "$wpep\t$hpep\n";
            }
        }
    }
}

