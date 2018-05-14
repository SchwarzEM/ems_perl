#!/usr/bin/env perl

# scaf2uniref_taxa.pl -- Erich Schwarz <emsch@caltech.edu>, 3/26/2012.
# Purpose: given a BlastP format-6 search of Uniref proteins with AUGUSTUS predictions, and the corresponding proteome, tabulate the taxa found for each scaffold.

use strict;
use warnings;
use Getopt::Long;

my $proteome = q{};
my $blastp   = q{};

my $scaffold = q{};
my $protein  = q{};
my $taxon    = q{};

my $data_ref;

my $help;

GetOptions ( 'proteome=s' => \$proteome,
             'blastp=s'   => \$blastp,
             'help'       => \$help     );

if ( $help or (! $proteome) or (! $blastp) ) { 
    die "Format: scaf2uniref_taxa.pl\n",
        "        --proteome|-p  [Uniref proteome used for BlastP]\n",
        "        --blastp|-b    [BlastP, in format-6 table, of AUGUSTUS proteome vs. Uniref]\n",
        "        --help|-h      [print this message]\n",
        ;
}

open my $PROTEOME, '<', $proteome or die "Can't open Uniref proteome $proteome: $!";
while (my $input = <$PROTEOME>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+) \s .+ Tax= (.+) [ ] RepID /xms ) { 
            $protein = $1;
            $taxon   = $2;
            $data_ref->{'protein'}->{$protein}->{'taxon'} = $taxon;
        }
        else { 
            die "Can't parse proteome entry: $input\n";
        }
    }
}
close $PROTEOME or die "Can't close filehandle to Uniref proteome $proteome: $!";

open my $BLASTP, '<', $blastp or die "Can't open BlastP $blastp $!";
while (my $input = <$BLASTP>) {
    chomp $input;  
    if ( $input =~ / /xms ) { 
        if ( $input =~ /\A (\S+)\.g\d+\.t\d+ \s+ (\S+) \s /xms ) {
            $scaffold = $1;
            $protein  = $2;
            $taxon    = $data_ref->{'protein'}->{$protein}->{'taxon'};
            $data_ref->{'scaffold'}->{$scaffold}->{'taxon'}->{$taxon}++;
        }
        else {
            die "Can't parse BlastP entry: $input\n";
        }
    }
}
close $BLASTP or die "Can't close filehandle to BlastP $blastp $!";

my @scaffolds = sort keys %{ $data_ref->{'scaffold'} };
foreach my $scaffold1 (@scaffolds) { 
    my @taxa = map { "$_ [$data_ref->{'scaffold'}->{$scaffold1}->{'taxon'}->{$_}]" } sort keys %{ $data_ref->{'scaffold'}->{$scaffold1}->{'taxon'} };
    my $tax_line = join '; ', @taxa;
    print "$scaffold1\t$tax_line\n";
}

