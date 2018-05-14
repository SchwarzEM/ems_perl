#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $srr    = q{};
my $prefix = q{};

my $header = '#!/bin/bash';
my $footer = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A .+ (SRR \d+) \s* \z/xms ) { 
        $srr = $1;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
    if ( $srr =~ /\A (SRR\d{3}) \d+ \z/xms ) { 
        $prefix = $1;
    }
    else {
        die "Can't parse putative SRR accession number : $srr\n";
    }
    print "$header\n\n" if $header;
    $header = q{};
    $footer = "\n";

                  # As described in "Determining the location of SRA data files for automated or scripted downloads."
                  #    in https://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.downloading_sra_data_using --
                  # 
                  # ftp://ftp-trace.ncbi.nih.gov
                  #                             /sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}
                  #     /<first 6 characters of accession>/<accession>/<accession>.sra

    print '    wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/' . $prefix . q{/} . $srr . q{/} . "$srr.sra ;\n";
}

print $footer;

