#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %sp2class = (
    'Acyrthosiphon pisum'  => 'Arthropods',
    'Apis mellifera'       => 'Arthropods',
    'Ceratitis capitata'   => 'Arthropods',
    'Diaphorina citri'     => 'Arthropods',
    'Hydra vulgaris'       => 'Metazoa',
    'Planococcus citri'    => 'Arthropods',
    'Trichoplax adhaerens' => 'Metazoa',
);

my %sp2prefix = (
    'Acyrthosiphon pisum'  => 'Arth_',
    'Apis mellifera'       => 'Arth_',
    'Ceratitis capitata'   => 'Arth_',
    'Diaphorina citri'     => 'Arth_',
    'Hydra vulgaris'       => 'Mzoa_',
    'Planococcus citri'    => 'Arth_',
    'Trichoplax adhaerens' => 'Mzoa_'
);  


while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ([A-Z][a-z]+ [_] (\S+)) .+ \[ (.+) \] \s* \z/xms ) {
            my $full_name = $1;
            my $orig_name = $2;
            my $species   = $3;
            my $class     = $sp2class{$species};
            my $sec_pref  = $sp2prefix{$species};
            print "$sec_pref$full_name\t$class\t$species\t$orig_name\t$orig_name\n";
        }
        else { 
            die "Cannot parse header: $input\n";
        }
    }
}


