#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my @seqs = ();

while (my $input = <>) { 
    chomp $input;
    # Sample input line:
    # ' CRE29430 (mss-1)   MIRKTTILFLAL-ALLATVARG-----ADGDNVEAGDAQLPSGGA--------------Q'
    if ( $input =~ /\A (.{20}) (\S+) \z/xms ) { 
        my $header   = $1;
        my $residues = $2;
        $header =~ s/\A\s+//;
        $header =~ s/\s+\z//;
        # Keep the order of sequences:
        push @seqs, $header; 
        $data_ref->{'seq'}->{$header}->{'residues'} .= $residues;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

foreach my $header (@seqs) {
    my $residues = $data_ref->{'seq'}->{$header}->{'residues'};
    print ">$header\n";
    my @output_lines 
        = unpack("a60" x (length($residues)/60 + 1), $residues);
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

