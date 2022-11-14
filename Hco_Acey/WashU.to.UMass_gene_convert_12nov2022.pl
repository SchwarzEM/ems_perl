#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw{ uniq };

my $washu2umass = q{};
my $genelist     = q{};

$washu2umass = $ARGV[0] if $ARGV[0];
$genelist    = $ARGV[1] if $ARGV[1];

my @umass_genes = ();

my $data_ref;

if ( (! -e $washu2umass ) or (! -e $genelist ) ) {
    die "Format: WashU.to.UMass_gene_convert_12nov2022.pl [WashU to UMass] [WashU gene list] > [mapped UMass gene list]\n";
}

open my $WASHU2UMASS, '<', $washu2umass;
while (my $input = <$WASHU2UMASS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $washu = $1;
        my $umass = $2;
        $data_ref->{'washu'}->{$washu}->{'umass'}->{$umass} = 1;
        $data_ref->{'umass'}->{$umass}->{'washu'}->{$washu} = 1;
    }
    else {
        die "In WashU to UMass $washu2umass, cannot parse input: $input\n";
    }
}
close $WASHU2UMASS;

open my $GENELIST, '<', $genelist;
while (my $input = <$GENELIST> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $gene = $1;
        my @umass_genes_new = ();
        if ( exists $data_ref->{'washu'}->{$gene}->{'umass'} ) {
            @umass_genes_new = sort keys %{ $data_ref->{'washu'}->{$gene}->{'umass'} };
            push @umass_genes, @umass_genes_new;
        }
    }
    else {
        die "In raw list $genelist, cannot parse input: $input\n";
    }
}
close $GENELIST;

@umass_genes = uniq(@umass_genes);
@umass_genes = sort(@umass_genes);

foreach my $umass_gene (@umass_genes) {
    print "$umass_gene\n";
}


