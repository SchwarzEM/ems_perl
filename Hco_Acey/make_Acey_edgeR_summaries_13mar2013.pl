#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

print '#!/bin/bash', "\n\n";
print "    rm Acey_edgeR_summaries_13mar2013.txt ;\n";
print "    touch Acey_edgeR_summaries_13mar2013.txt ;\n\n";
print "    echo >> Acey_edgeR_summaries_13mar2013.txt ;\n";

while (my $input = <>) {
    chomp $input;
    my $file_base = basename $input;
    my $label     = $file_base;

    if ( $label =~ /\.DGEList\.exactTest\.topTags\.txt/ ) { 
        die "Looks like bad input\n";
    }

    $label =~ s/\.edgeR\.txt//;
    $label =~ s/_and_/ and /g;
    $label =~ s/\.vs\./ vs /g;

    print "    echo \"$label\" >> Acey_edgeR_summaries_13mar2013.txt ;\n";
    print "    echo \"Upregulated, q <= 1e-03\" >> Acey_edgeR_summaries_13mar2013.txt ;\n";
    print "    cat $input ",
          ' | grep Acey_ | perl -ne \' $input = $_ ; chomp $input; if ( ( $input =~ /\A\S+\s+(\S+)\s+(\S+)\z/ ) and ( $2 <= 1e-03 ) and ( $1 > 0 ) ) { print "$input\n"; } \' | wc -l',
          " >> Acey_edgeR_summaries_13mar2013.txt ;\n",
          ;
    print "    echo \"Downregulated, q <= 1e-03\" >> Acey_edgeR_summaries_13mar2013.txt ;\n"; 
    print "    cat $input ",
          ' | grep Acey_ | perl -ne \' $input = $_ ; chomp $input; if ( ( $input =~ /\A\S+\s+(\S+)\s+(\S+)\z/ ) and ( $2 <= 1e-03 ) and ( $1 < 0 ) ) { print "$input\n"; } \' | wc -l',
          " >> Acey_edgeR_summaries_13mar2013.txt ;\n",
          ;
    print "    echo >> Acey_edgeR_summaries_13mar2013.txt ;\n";
    print "\n";
}

