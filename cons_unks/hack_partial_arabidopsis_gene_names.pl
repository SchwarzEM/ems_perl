#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref; 

my %part2full   = ();
my @input_lines = ();

# Note that this does not even try to keep track of whether there are promiscuous links of TAIR ID to human-readble name;
#    the last linkage recorded for a given 'part' is what stands.  This is a hack.

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \t (\S+) \b .* \z/xms ) { 
        my $gene_name = $1;
        if ( $gene_name =~ /[|]/xms ) {
            my @parts = split /[|]/, $gene_name;
            foreach my $part (@parts) { 
                $part2full{$part} = $gene_name;
            }
        }
    }
    else { 
        die "Can't parse input line: $input\n";
    }
    push @input_lines, $input;
}

foreach my $output (@input_lines) { 
    if ( $output =~ /\A (\S+) \t (\S+) (\b .*) \z/xms ) {
        my $uniprot   = $1;
        my $gene_name = $2;
        my $desc      = $3;
        if ($part2full{$gene_name}) {
            $gene_name = $part2full{$gene_name};
        }
        print "$uniprot\t$gene_name$desc\n";
    }
    else { 
        die "Can't parse output line: $output\n";
    }
}

