#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my $gene = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A Sequence \s+ \" (\S+\.g\d+) \.t\d+ \" /xms ) { 
        $gene = $1;
    }
    elsif ( $input =~ /\A InterPro \s+ (IPR\d+) \s+ (\S.+\S) \s* \z/xms ) { 
        my $desc = $2;
        my $ipr  = $1;
        my $full_desc = "$desc [$ipr]";
        $data_ref->{'gene'}->{$gene}->{'interpro'}->{$full_desc} = 1;
    }
}

my $header = "Gene\tInterPro\n";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @ipr_descs     = ();
    my $ipr_desc_text = q{};
    if ( exists $data_ref->{'gene'}->{$gene1}->{'interpro'} ) { 
        @ipr_descs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'interpro'} };
        $ipr_desc_text = join '; ', @ipr_descs;
    }
    print $header if $header;
    $header = q{};
    print "$gene1\t$ipr_desc_text\n";
}

