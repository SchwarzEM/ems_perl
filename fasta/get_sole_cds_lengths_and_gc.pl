#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $seq = q{};
my $data_ref;

my $header = "\tLength\tGC";

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+) /xms ) { 
            $seq = $1;
            if ( exists $data_ref->{'seq'}->{$seq} ) { 
                die "Redundant sequence name: $seq\n";
            }
        }
        else { 
            die "Can't parse input line: $input\n";
        }
    }
    elsif ( $input =~ /\S/xms ) { 
            if (! $seq) { 
                die "Cannot map to a named sequence: $input\n";
            }
            $input =~ s/\s//g;
            $data_ref->{'seq'}->{$seq}->{'residues'} .= $input;
    }
}

my @seqs = sort keys %{ $data_ref->{'seq'} };
if (@seqs) {
    print "$header\n";
}
foreach my $seq1 (@seqs) { 
    my $residues = $data_ref->{'seq'}->{$seq1}->{'residues'};
    my $length   = length($residues);
    my $gc_count = ($residues =~ tr/cCgG/cCgG/);
    my $at_count = ($residues =~ tr/aAtT/aAtT/);
    my $nt_count = $gc_count + $at_count;
    my $gc_frac  = ( $gc_count / $nt_count );
    $gc_frac     = sprintf("%.3f", $gc_frac);
    print "$seq1\t$length\t$gc_frac\n";
}


