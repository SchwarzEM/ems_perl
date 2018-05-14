#!/usr/bin/env perl

use strict;
use warnings;

my $gene2go = 'gene_association.WS230.wb.ce';

if (! -e $gene2go) { 
    die "Can't find gene2go $gene2go\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+)\.txt \z/xms ) { 
        my $stem   = $1;
        my $output = $stem . '.func_input.txt';
        $output    = safename($output);
        print "    qual2func_table.pl -g $gene2go -t $input > $output ;\n";
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

