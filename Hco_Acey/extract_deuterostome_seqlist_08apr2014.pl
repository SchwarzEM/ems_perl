#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %prefix2spp = (
    'Bflo' => 'Branchiostoma floridae',
    'Cfam' => 'Canis familiaris',
    'Cint' => 'Ciona intestinalis',
    'Dnov' => 'Dasypus novemcinctus',
    'Drer' => 'Danio rerio',
    'Eeur' => 'Erinaceus europaeus',
    'Ggal' => 'Gallus gallus',
    'Hsap' => 'Homo sapiens',
    'Lafr' => 'Loxodonta africana',
    'Lcha' => 'Latimeria chalumnae',
    'Maur' => 'Mesocricetus auratus',
    'Mdom' => 'Monodelphis domestica',
    'Mmus' => 'Mus musculus',
    'Oana' => 'Ornithorhynchus anatinus',
    'Pmar' => 'Petromyzon marinus',
    'Skow' => 'Saccoglossus kowalevskii',
    'Spur' => 'Strongylocentrotus purpuratus',
    'Tgut' => 'Taeniopygia guttata',
    'Xtro' => 'Xenopus tropicalis',
);

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > ( ([A-Z][a-z]{3}) _ \S+) /xms ) {
        my $seqname = $1;
        my $prefix  = $2;
        if ( exists $prefix2spp{$prefix} ) {
            print "$seqname\n";
        }
    }
}

