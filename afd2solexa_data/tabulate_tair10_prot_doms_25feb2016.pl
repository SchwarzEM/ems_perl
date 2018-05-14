#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tDomains";

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+) \.\d+ \t (?: [^\t]* \t){9} (IPR\d+) \t ([^\t]+) \z/xms ) { 
        my $gene     = $1;
        my $ipr_id   = $2;
        my $ipr_desc = $3;
        if ( $ipr_desc =~ /\A ([^;]+) [;] /xms ) { 
            $ipr_desc = $1;
        }
        $ipr_desc  =~ s/\A\s+//;
        $ipr_desc  =~ s/\s+\z//;
        my $domain = "$ipr_desc [$ipr_id]";
        $data_ref->{'gene'}->{$gene}->{'domain'}->{$domain} = 1;
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @domains = sort keys %{ $data_ref->{'gene'}->{$gene}->{'domain'} };
    my $domain_text = join '; ', @domains;
    print "$header\n" if $header;
    $header = q{};
    print "$gene\t$domain_text\n";
}

