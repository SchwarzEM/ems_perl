#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = $ARGV[0];

my %taxa      = ();
my %multigene = ();

# Read the file once to see which species need to get serial numbers of genes;
#     then, reread it to make those gene numbers happen.

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A > .+ \[ (\S+) \s+ (\S+) .* \] \s* \z/xms ) {
        my $genus   = $1;
        my $species = $2;
        my $taxon   = $genus. q{_} . $species;
        $taxa{$taxon}++;
        my $gene_no =  $taxa{$taxon};
        if ( $gene_no >= 2 ) {
            $multigene{$taxon} = 1;
        }
    }
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse: $input\n";
    }
}
close $INFILE;

# Result the gene-counter to zero:
%taxa = ();

open $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A > ( .+ \[ (\S+) \s+ (\S+) .* \] ) \s* \z/xms ) {
        my $header  = $1;
        my $genus   = $2;
        my $species = $3;
        my $taxon   = $genus. q{_} . $species;
        $taxa{$taxon}++;
        my $gene_no =  $taxa{$taxon};
        my $name = q{};

        if ( $multigene{$taxon} ) {
            $name = "$taxon.$gene_no";
        }
        else {
            $name = $taxon;
        }
        print ">$name\t$header\n";
    }
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse: $input\n";
    }
    else {
        print "$input\n";
    }
}
close $INFILE;

