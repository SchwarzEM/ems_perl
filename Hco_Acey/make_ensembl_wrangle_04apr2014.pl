#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

while (my $input = <>) { 
    # Sample input:
    # ../seqs/Canis_familiaris.CanFam3.1.75.pep.all.fa
    chomp $input;
    my $basename = basename($input);
    if ( $basename =~ /\A (([A-Z]) [a-z]+ _ ([a-z]{3}) [a-z]+ \. \S+) \.pep\.all\.fa \z/xms ) { 
        my $stem    = $1;
        my $sp_init = $2;
        my $sp_3let = $3;
        my $prefix  = $sp_init . $sp_3let . q{_};
        my $output  = 'proteomes/' . $stem . '.nice.pep.all.fa';
        print "    get_largest_isoforms.pl -i $input -t ens | prettify_ensembl_proteome_headers_18apr2014.pl $prefix - > $output ;\n";
    }
    else { 
        die "Can't parse input filename: $input\n";
    }
}

