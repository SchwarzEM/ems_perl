#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $allele_list = $ARGV[0];
my $pheno_table = $ARGV[1];

my $data_ref;

open my $ALLELE_LIST, '<', $allele_list;
while (my $input = <$ALLELE_LIST> ) {
    chomp $input;
    $data_ref->{'banned_allele'}->{$input} = 1;
}
close $ALLELE_LIST;

open my $PHENO_TABLE, '<', $pheno_table;
while (my $input = <$PHENO_TABLE> ) {
    chomp $input;
    if ( $input =~ /\A (?: [^\t]* \t){5} (\S+) \t/xms ) { 
        my $allele = $1;
        if (! exists $data_ref->{'banned_allele'}->{$allele} ) { 
            print "$input\n";
        }
    }
}
close $PHENO_TABLE;


