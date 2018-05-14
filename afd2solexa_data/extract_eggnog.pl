#!/usr/bin/env perl

# extract_eggnog.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/23/2010.
# Purpose: extract compact C. elegans-specific annotations, one line per gene (sequence name only), from eggNOG 2.0 orthgroups.mapping.txt.

use strict;
use warnings;

my $data_ref;

my $cds    = q{};
my $kog   = q{};
my $annot = q{};

# 6239.C17G10.9a.2        32      537     meNOG04182      Eukaryotic translation initiation factor 3, subunit
# 6239.C17G10.9a.2        1       537     KOG3677 RNA polymerase I-associated factor - PAF67

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A 6239 \. ([^\.]+ \. \d+) \S* \t \d+ \t \d+ \t (\S+OG\S+) \t (.+) \z/xms ) { 
        $cds   = $1;
        $kog   = $2;
        $annot = $3;
        my $text = "$kog: $annot";       
        $data_ref->{'gene'}->{$cds}->{'kog'}->{$text} = 1;
    }
    elsif ( $input =~ / \A 6239 \. /xms ) { 
        die "Can't parse: $input\n";
    }
}

foreach my $g1  (sort keys %{ $data_ref->{'gene'} } ) { 
    my @annots = sort keys %{ $data_ref->{'gene'}->{$g1}->{'kog'} };
    my $annot_line = join '; ', @annots;
    print "$g1\t\"$annot_line\"\n";
}

