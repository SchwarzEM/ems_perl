#!/usr/bin/env perl

use strict;
use warnings;

if ( (! -e 'LSC_run_28nov2012/corrected' ) or (! -e 'LSC_run_28nov2012/uncorrected' ) or (! -e 'LSC_run_28nov2012/full' ) ) {
    die "Can't find destination directory for renamed files\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A LSC_run(\d+)\/output\.\d+\/(\S+)_LR_SR\.map\.fa \z/xms ) { 
        my $index    = $1;
        my $type     = $2;
        my $new_name = 'Csp9_LSC_run_2012.11.28.' . $index . q{_} . $type . '.fa';
        print "    mv -i $input LSC_run_28nov2012/$type/$new_name ;\n";
    }
    else { 
        die "Can't parse: $input\n";
    }
}

