#!/usr/bin/env perl

use strict;
use warnings;

my $gene2go = '../../augustus_preds_24oct2012/blast2go/Acey_2012.10.24.pep.max_isos.blast2go.annot.filt.txt';

if (! -e $gene2go) { 
    die "Can't find gene2go $gene2go\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (Acey_RSEM_01dec2012\.\S+)\.txt \z/xms ) { 
        my $stem   = $1;
        my $output = $stem . '.func_input.txt';
        $output    = safename($output);
        print "    qual2func_table_25nov2012.pl -g $gene2go -t $input > $output ;\n";
    }
    else { 
        die "Can't parse input $input\n";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

