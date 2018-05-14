#!/usr/bin/perl

# replace_clone_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/2/2007.
# Purpose: Go through a Kuntz-related text document and systematically replace old with new clone names.

use strict;
use warnings;

my %renamed = ( 
    'CB5161_anon1.tfa'         => "Cbre_JD01",
    'CB5161_anon2.tfa'         => "Cbre_JD02",
    'CB5161_anon3A.tfa'        => "Cbre_JD03",
    'CB5161_anon3B.tfa'        => "Cbre_JD04",
    'CB5161_anon4A.tfa'        => "Cbre_JD05",
    'CB5161_anon4B.tfa'        => "Cbre_JD06",
    'CB5161_ceh-13.lin-39.tfa' => "Cbre_JD07",
    'CB5161_ceh-8.tfa'         => "Cbre_JD08",
    'CB5161_col-14.tfa'        => "Cbre_JD09",
    'CB5161_egl-44.tfa'        => "Cbre_JD10",
    'CB5161_egl-46.tfa'        => "Cbre_JD11",
    'CB5161_gpa-6.tfa'         => "Cbre_JD12",
    'CB5161_lin-11.tfa'        => "Cbre_JD13",
    'CB5161_lin-3.tfa'         => "Cbre_JD14",
    'CB5161_ndescA.tfa'        => "Cbre_JD15",
    'CB5161_ndescB.tfa'        => "Cbre_JD16",
    'CB5161_nlp-8.tfa'         => "Cbre_JD17",
    'CB5161_rde-1A.tfa'        => "Cbre_JD18",
    'CB5161_rde-1B.tfa'        => "Cbre_JD19",
    'CB5161_ref-1.tfa'         => "Cbre_JD20",
    'CB5161_srw-2.tfa'         => "Cbre_JD21",
    'CB5161_weird.tfa'         => "Cbre_JD22",

    'PS1010_anonymous.tfa'     => "Csp3_JD01",
    'PS1010_ceh-13.lin-39.tfa' => "Csp3_JD02",
    'PS1010_egl-30.tfa'        => "Csp3_JD03",
    'PS1010_egl-5.mab-5.tfa'   => "Csp3_JD04",
    'PS1010_lin-11.tfa'        => "Csp3_JD05",
    'PS1010_lin-3.tfa'         => "Csp3_JD06",
    'PS1010_php-3.tfa'         => "Csp3_JD07",
);

my @old_names = keys %renamed;

while (my $input = <>) { 
    foreach my $old_name (@old_names) { 
        if ( $input =~ /$old_name/xms ) { 
            $input =~ s/$old_name/$renamed{$old_name}/g;
        }
    }
    print $input;
}

