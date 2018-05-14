#!/usr/bin/env perl

use strict;
use warnings;

my $header  = "Gene\tChr_coords";
my @outputs = ();

while (my $input = <>) {
    chomp $input;

    # Sample input:
    # Gene WB ID	Gene Public Name	Chr Name	Strand	Start (bp)	End (bp)
    # WBGene00000001	aap-1	I	1	5107844	5110167

    if ($input =~ /\A (WBGene\d+) \t [^\t]* \t (\S+) \t (\S+) \t (\d+) \t (\d+) /xms ) { 
        my $gene    = $1;
        my $chr     = $2;
        my $ori_val = $3;
        my $start   = $4;
        my $stop    = $5;
        my $ori     = q{};
        if ( $ori_val == 1 ) {
            $ori = q{+};
        }
        elsif ( $ori_val == -1 ) { 
            $ori = q{-};
        }
        my $out_text = $gene . "\t" . $chr . q{:} . $start . q{-} . $stop . q{ [} . $ori . q{]};
        push @outputs, $out_text; 
    }
}

foreach my $output (@outputs) { 
    print "$header\n" if $header;
    $header = q{};
    print "$output\n";
}

