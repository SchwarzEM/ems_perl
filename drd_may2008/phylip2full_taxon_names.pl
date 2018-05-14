#!/usr/bin/env perl

# phylip2full_taxon_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/13/2008.
# Purpose: fix stupid munging of >10-character taxon names in a PHYLIP tree, given the source FASTA.

use strict;
use warnings;
use Getopt::Long;

my $fasta          = q{};
my $phylip_tree    = q{};
my %phyl2full_name = ();

GetOptions ( "fasta=s"  => \$fasta,  
             "phylip=s" => \$phylip_tree, );

if ( (! -r $fasta) or (! -r $phylip_tree ) ) { 
    die "Format: ./phylip2full_taxon_names.pl",
        " --fasta=[FASTA w/ full names]",
        " --phylip=[PHYLIP tree with munged names]",
        "\n",
        ;
}

open my $FASTA, '<', $fasta or die "Can't open FASTA $fasta: $!";

while (my $input = <$FASTA>) { 
    if ( $input =~ /\A > ( \S+ ) /xms ) { 
        my $full_name = $1;
        my $phyl_name = substr($full_name, 0, 10);
        $phyl2full_name{$phyl_name} = $full_name;
    }
}
close $FASTA or die "Can't close filehandle for FASTA $fasta: $!";

my @phylip_names = sort keys %phyl2full_name;

open my $PHYLIP, '<', $phylip_tree 
    or die "Can't open PHYLIP tree $phylip_tree: $!";

# Computationally inefficient loop, but should at least work reliably:
while (my $input = <$PHYLIP>) { 
    foreach my $phyl_name (@phylip_names) { 
        my $full_name = $phyl2full_name{$phyl_name};
        $input =~ s/$phyl_name/$full_name/g;
    }
    print $input;
}

close $PHYLIP or die "Can't close filehandle for PHYLIP tree $phylip_tree: $!";

