#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $list   = q{};
my $spans  = q{};
my %select = ();

$list  = $ARGV[0] if $ARGV[0];
$spans = $ARGV[1] if $ARGV[1];

if ( (! $list ) or (! $spans ) ) {
    die "Format: filter_aug_gene_lines_24jun2018.pl [gene list] [AUGUSTUS-GFF all-gene-lines file] > [only GFF lines in list] ;\n";
}

open my $LIST, '<', $list;
while (my $gene = <$LIST>) { 
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        $select{$gene} = 1;
    }
    else {
        die "From gene-list file $list, cannot parse: $gene\n";
    }
}
close $LIST;

# Sample input line:
# chrI_pilon      AUGUSTUS        gene    6054    11981   0.98    -       .       ID=chrI_pilon.g1

open my $SPANS, '<', $spans;
while (my $input = <$SPANS>) {
    chomp $input;
    if ( $input =~ /\A .+ \t ID=(\S+) \b [^\t]* \z /xms ) {
        my $gene = $1;
        if ( exists $select{$gene} ) {
            print "$input\n";
        }
    }
    else {
        print "From gene-spans file $spans, cannot parse: $input\n";
    }
}
close $SPANS;

