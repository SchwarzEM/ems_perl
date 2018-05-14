#!/usr/bin/env perl

use strict;
use warnings;

my %synonyms = ( 
    'Acey_2012.08.05_0248.g90' => 'ASP-1',
    'Acey_2012.08.05_0184.g991' => 'ASP-2',
    'Acey_2012.08.05_0067.g52' => 'ASP-3A',
    'Acey_2012.08.05_0067.g57' => 'ASP-3B',
    'Acey_2012.08.05_0004.g1750' => 'ASP-4A',
    'Acey_2012.08.05_0004.g1747' => 'ASP-4B',
    'Acey_2012.08.05_0003.g1624' => 'ASP-5A',
    'Acey_2012.08.05_0003.g1625' => 'ASP-5B',
    'Acey_2012.08.05_0067.g36' => 'ASP-6',
    'Acey_2012.08.05_0222.g2615' => 'ASP-7A',
    'Acey_2012.08.05_0012.g1631' => 'ASP-7B',
    'Acey_2012.08.05_0067.g64' => 'ASP-NIF-A',
    'Acey_2012.08.05_0067.g62' => 'ASP-NIF-B',
);

# >ASP-s0001.g198
# >ASPR-s0002.g551

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > (Acey_2012\.08\.05_(\d+\.g\d+))\.t\d+ /xms ) { 
            my $cds   = $1;
            my $label = $2;
            my $gene  = 'ASP-s' . $label;
            if ( exists $synonyms{$cds} ) {
                $gene = $synonyms{$cds};
            }
            print ">$gene\n";
        }
        else {
            die "Can't parse header: $input\n";
       }
    }
    else { 
        print "$input\n";
    }
}

