#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $category = $ARGV[0];
my $tax_ids  = $ARGV[1];
my $headers  = $ARGV[2];

my %ok_taxid = ();

if ( $category !~ /\A [A-Za-z]+ \z/xms ) { 
   die "Do you really want this as a category tag?: $category\n";
}

open my $TAX_IDS, '<', $tax_ids;
while (my $input = <$TAX_IDS>) {
    chomp $input;
    if ( $input !~ /\A \d+ \z/xms ) { 
        die "Can't parse alleged tax ID number from file $tax_ids: $input\n";
    }
    $ok_taxid{$input} = 1;
}
close $TAX_IDS;

open my $HEADERS, '<', $headers;
while (my $input = <$HEADERS>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) \b .* \s NCBI_TaxID = (\d+) \s /xms ) {
        my $protein = $1;
        my $taxid   = $2;
        if ( exists $ok_taxid{$taxid} ) { 
            print "$category\t$protein\n";
        }
    }
}
close $HEADERS;


