#!/usr/bin/env perl

use strict;
use warnings;

if (! -e 'LSC_run_28nov2012/corrected') { 
    die "Can't find target directory.\n";
}

print "\n";
print "    rm LSC_run_28nov2012/Csp9_LSC_run_2012.11.28_corrected.fa ;\n";
print "    touch LSC_run_28nov2012/Csp9_LSC_run_2012.11.28_corrected.fa ;\n";
print "\n";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (LSC_run_28nov2012\/corrected\.orig\/(Csp9_LSC_run_2012\.11\.28\.\d+))_corrected\.fa \z/xms ) { 
        my $orig_file_stem = $1;
        my $orig_name_stem = $2;
        my $orig_file      = $orig_file_stem . '_corrected.fa';
        my $orig_rename    = $orig_file_stem . '_corrected.orig.fa';
        print "    mv -i $orig_file $orig_rename ;\n";
        print "    tag_FASTA_names.pl -p $orig_name_stem", "_ $orig_rename >> LSC_run_28nov2012/Csp9_LSC_run_2012.11.28_corrected.fa ;\n", ;
        print "\n";
    }
    else { 
        die "Can't parse: $input\n";
    }
}

