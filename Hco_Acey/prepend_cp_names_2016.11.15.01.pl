#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %cds2gene = ( 
    'Acey_s0154.g3007.t1' => 'CP-1',
    'Acey_s0154.g3016.t2' => 'CP-2',
    'Acey_s0154.g3018.t8' => 'CP-3',
    'Acey_s0154.g3005.t2' => 'CP-4',
    'Acey_s0619.g723.t1'  => 'CP-5',
    'Acey_s0571.g111.t1'  => 'CP-6',
    'Acey_s0154.g2963.t1' => 'CP-7',
    'Acey_s0028.g1832.t1' => 'CP-8',
    'Acey_s0220.g2503.t1' => 'CP-9',
    'Acey_s0007.g3469.t1' => 'CP-10',
    'Acey_s0040.g229.t3'  => 'CP-11',
    'Acey_s0619.g718.t2'  => 'CP-12',
);

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A [>](\S+) (.*)\z/xms ) { 
        my $cds   = $1;
        my $annot = $2;
        my $gene  = $cds2gene{$cds};
        $input = '>' . $gene . q{_} . $cds . $annot;
    }
    print "$input\n";
}


