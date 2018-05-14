#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (.*) \z/xms ) {  
        my $uniprot = $1;
        my $gene    = $2;
        my $desc    = $3;
        my $alias   = q{};
        if ( $gene =~ /\A ENSG\d+ \| (\S+) \z/xms ) {
            $alias = $1;
            $data_ref->{'alias'}->{$alias}->{'gene'}              = $gene;
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
            if ( $desc =~ /\S/xms ) { 
                $desc =~ s/["]//g;
                $data_ref->{'gene'}->{$gene}->{'desc'}->{$desc} = 1;
            }
        }
        else { 
            $data_ref->{'odd_uniprot'}->{$uniprot}->{'odd_gene'}->{$gene} = 1;
            # Map this *both* ways, so that carryover is possible later on.
            $data_ref->{'odd_gene'}->{$gene}->{'odd_uniprot'}->{$uniprot} = 1;
            if ( $desc =~ /\S/xms ) {
                $desc =~ s/["]//g;
                $data_ref->{'odd_gene'}->{$gene}->{'desc'}->{$desc} = 1;
            }
        }
    }
    else { 
        die "Can't parse input: $input\n";
    }
}

# If possible, carry over any proteins and descriptions from odd to normal genes.
my @odd_genes = sort keys %{ $data_ref->{'odd_gene'} };
foreach my $odd_gene (@odd_genes) {
    if ( exists $data_ref->{'alias'}->{$odd_gene}->{'gene'} ) {
        my $true_name = $data_ref->{'alias'}->{$odd_gene}->{'gene'};

        # Uniprot carryover:
        if ( exists $data_ref->{'odd_gene'}->{$odd_gene}->{'odd_uniprot'} ) { 
            my @odd_uniprots = sort keys %{ $data_ref->{'odd_gene'}->{$odd_gene}->{'odd_uniprot'} };
            foreach my $odd_uniprot (@odd_uniprots) { 
                $data_ref->{'uniprot'}->{$odd_uniprot}->{'gene'}->{$true_name} = 1;                
                delete $data_ref->{'odd_gene'}->{$odd_gene}->{'odd_uniprot'}->{$odd_uniprot};
                delete $data_ref->{'odd_uniprot'}->{$odd_uniprot}->{'odd_gene'}->{$odd_gene};
            }
        }

        # Description carryover:
        if ( exists $data_ref->{'odd_gene'}->{$odd_gene}->{'desc'} ) {
             my @odd_gene_descs = sort keys %{ $data_ref->{'odd_gene'}->{$odd_gene}->{'desc'} };
             foreach my $odd_gene_desc (@odd_gene_descs) { 
                 $data_ref->{'gene'}->{$true_name}->{'desc'}->{$odd_gene_desc} = 1;
             }
        }
    }
}

my @first_uniprots = sort keys %{ $data_ref->{'uniprot'} };
foreach my $uniprot (@first_uniprots) { 
    my @first_genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene'} };
    foreach my $gene (@first_genes) {
        my $desc = export_desc_set('gene', $gene);
        print "$uniprot\t$gene\t$desc\n";
    }
}

my @odd_uniprots = sort keys %{ $data_ref->{'odd_uniprot'} };
foreach my $uniprot (@odd_uniprots) {
    my @odd_genes_local = sort keys %{ $data_ref->{'odd_uniprot'}->{$uniprot}->{'odd_gene'} };
    foreach my $odd_gene_local (@odd_genes_local) { 
        my $desc = export_desc_set('odd_gene', $odd_gene_local);
        print "$uniprot\t$odd_gene_local\t$desc\n";
    }
}

sub export_desc_set {
    my $_gene_type = $_[0];
    my $_gene_name = $_[1];
    my $_desc      = q{};
    if ( exists $data_ref->{$_gene_type}->{$_gene_name}->{'desc'} ) {
        my @_descs = sort keys %{ $data_ref->{$_gene_type}->{$_gene_name}->{'desc'} };
        $_desc = join '; ', @_descs;
        $_desc = "\"$_desc\"";
    }
    return $_desc;
}
