#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $es_list = q{};
my $uni2par = q{};

$es_list = $ARGV[0] if $ARGV[0];
$uni2par = $ARGV[1] if $ARGV[1];

my %uni2para = ();

if ( (! -e $es_list ) or (! -e $uni2par ) ) {
    die "Format: angio_es_map_05feb2023.pl [ES UniProt list] [UniProt to ParaSite table] > [ES ParaSite gene list];\n";
}

open my $UNI2PAR, '<', $uni2par;
while ( my $input = <$UNI2PAR>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t \S+ \z/xms ) {
        my $uniprot       = $1;
        my $parasite_gene = $2;
        $uni2para{$uniprot} = $parasite_gene;
    }
    else {
        die "From UniProt to ParaSite table $uni2par, cannot parse: $input\n";
    }
}
close $UNI2PAR;

open my $ES_LIST, '<', $es_list;
while ( my $input = <$ES_LIST>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
       	my $es_uniprot = $1;
        if ( exists $uni2para{$es_uniprot} ) {
            print "$uni2para{$es_uniprot}\n";
        }
    }
    else {
       	die "From ES list $es_list, cannot parse: $input\n";
    }
}
close $ES_LIST;
