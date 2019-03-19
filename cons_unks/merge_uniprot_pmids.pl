#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]*) \t ([^\t]*) \z/xms ) {
        my $uniprot    = $1;
        my $pmid_text1 = $2;
        my $pmid_text2 = $3;
        my @pmids      = split /; /, $pmid_text1;
       	my @more_pmids = split /; /, $pmid_text2;
        push @pmids, @more_pmids;
        @pmids = sort @pmids;
        @pmids = uniq(@pmids);
        foreach my $pmid (@pmids) {
            if ( $pmid !~/ \A \d+ \z/xms ) {
                die "Invalid PMID ($pmid) in: $input\n";
            }
            $data_ref->{'uniprot'}->{$uniprot}->{'pmid'}->{$pmid} = 1;
        }
    }
    else {
        die "Cannot parse: $input\n";
    }
}

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
foreach my $uniprot (@uniprots) {
    my @pmids = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'pmid'} };
    my $pmid_text = join '; ', @pmids;
    print "$uniprot\t$pmid_text\n";
}

