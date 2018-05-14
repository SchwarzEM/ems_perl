#!/usr/bin/perl

use strict;
use warnings;

my $seqname = q{};
my $caps    = 0;
my $j       = 0;
my %names2residues_ref = ();

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A > (\S+) /xms) { 
        $seqname = $1;
    }
    else { 
        $input =~ tr/[^a-zA-Z]//;
        my @chars = split //, $input;
        push @{ $names2residues_ref{ $seqname } }, @chars;
    }
}

foreach my $name (sort keys %names2residues_ref) { 
    print "$name: ";
    my $length = @{ $names2residues_ref{ $seqname } };
    print "$length residues\n";
    foreach my $i (1..$length) { 
        $j = $i - 1;
        if ( ${ $names2residues_ref{ $seqname } }[$j] =~ /[A-Z]/ ) { 
            if (! $caps) { 
                print "$i-";
                $caps = $i;
            }
        }
        if ( ${ $names2residues_ref{ $seqname } }[$j] =~ /[a-z]/ ) {
            if ($caps) {
                print "$j\n";  # The *previous* residue was CAPs.
                $caps = 0;
            }
        }
    }
}

# Edge case: if sequence ends with a CAPs residue
if ($caps) {
    $j++;          # The *current* residue is CAPS!
    print "$j\n";
}

