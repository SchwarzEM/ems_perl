#!/usr/bin/env perl

# wbg_list2tab_annot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/17/2010.
# Purpose: given a gene list and a user-specified annotation, create a simple table of genes/annots.

use strict;
use warnings;
use Getopt::Long;

my $annot = q{};
my $help;

GetOptions ( 'annot=s', => \$annot,  
             'help'     => \$help, );

if ($help or (! $annot)) { 
    die "Format: wbg_list2tab_annot.pl --annot|-a [annotation text] <file or stream>\n";
}

if ( $annot =~ /\t/xms ) { 
    die "The annotation \"$annot\" contains a tab, which would split the data fields.\n";
}

while (my $input = <>) { 
    chomp $input;
    print "$input\t$annot\n";
}


