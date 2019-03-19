#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

# Sample input line:
# A0A068FL09      211     526     211     527     PF01501 Glyco_transf_8 [etc.]

while (my $input = <>) {
    chomp $input;
    if ( ( $input !~ /\A [#]/xms ) and ( $input =~ /\A (\S+) \t (?: [^\t]* \t){4} (\S+) \t/xms ) ) {
        my $uniprot = $1;
        my $pfam = $2;
        $data_ref->{'uniprot'}->{$uniprot}->{'pfam'}->{$pfam} = 1;
    }
    elsif ( $input !~ /\A [#]/xms ) {
        die "Cannot parse input: $input\n";
    }
}

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
foreach my $uniprot (@uniprots) {
    my @pfams = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'pfam'} };
    my $pfam_text = join ';', @pfams;
    print "$uniprot\t$pfam_text\n";
}

