#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $oma_uniprot  = q{};
my $uniprot2gene = q{};
my $oma_groups   = q{};

my $data_ref;

$oma_uniprot  = $ARGV[0] if $ARGV[0];
$uniprot2gene = $ARGV[1] if $ARGV[1];
$oma_groups   = $ARGV[2] if $ARGV[2];

if ( (! $oma_uniprot) or (! $uniprot2gene) or (! $oma_groups) ) {
     die "Format: uniprot2genes.oma_12nov2018.pl [oma-uniprot] [uniprot2gene] [oma-groups] > [gene2oma_groups]\n";
}

# Sample oma-uniprot.txt file:
# ARATH00001    NAC1_ARATH
# ARATH00001    Q0WV96
# ARATH00001    A0A178WAE4

open my $OMA_UNIPROT, '<', $oma_uniprot;
while (my $input = <$OMA_UNIPROT>) {
    chomp $input;
    if ( $input	=~ /\A (\S+) \t	(\S+) \z/xms ) { 
       	my $oma_id  = $1;
       	my $uniprot = $2;
       	$data_ref->{'oma_id'}->{$oma_id}->{'uniprot'}->{$uniprot} = 1;
    }
    else {
       	die "In	oma_uniprot file $oma_uniprot, cannot parse: $input\n";
    }
}
close $OMA_UNIPROT;

# Sample uniprot2gene.txt file:
# Q06469        yeast|S000006362|CUR1
# P25044        yeast|S000002389|PTP1
# P53298        yeast|S000003411|OKP1
# Q9FEF7	arabidopsis|AT1G24822; arabidopsis|AT1G25097

open my $UNIPROT2GENE, '<', $uniprot2gene;
while (my $input = <$UNIPROT2GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S[^\t]+\S) \z/xms ) { 
        my $uniprot   = $1;
        my $gene_text = $2;
        my @genes = split /[;] /, $gene_text;
        foreach my $gene (@genes) {
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
        }
    }
    else {
        die "In uniprot2gene file $uniprot2gene, cannot parse: $input\n";
    }
}
close $UNIPROT2GENE;

# Sample oma-groups.txt file:
# 183871        YREWNLG SCHPO04143	YEAST00244
# 188642        VWVIRRV SCHPO02325	YEAST03832
# 189800        AEWAYVP SCHPO05002	YEAST02520

open my $OMA_GROUPS, '<', $oma_groups;
while (my $input = <$OMA_GROUPS>) {
    chomp $input;

    if ( $input =~ /\A (\d+) \t (\S+) \t (.+) \z/xms ) {
        my $oma_no       = $1;
        my $oma_print    = $2;
        my $oma_id_text  = $3;

        my $oma_family   = "oma\|$oma_no\|$oma_print";

        my @oma_ids      = split /\t/, $oma_id_text;
        my @oma_uniprots = ();
        foreach my $oma_id (@oma_ids) {
            if ( exists $data_ref->{'oma_id'}->{$oma_id}->{'uniprot'} ) {
                my @new_oma_uniprots = sort keys %{ $data_ref->{'oma_id'}->{$oma_id}->{'uniprot'} };
                push @oma_uniprots, @new_oma_uniprots;
            }
            else {
                warn "Cannot map OMA ID $oma_id to UniProt ID in: $input\n";
            }
        }
        @oma_uniprots = sort @oma_uniprots;
        @oma_uniprots = uniq(@oma_uniprots);

        my @oma_genes = ();
        foreach my $oma_uniprot (@oma_uniprots) {
            if ( exists $data_ref->{'uniprot'}->{$oma_uniprot}->{'gene'} ) {
                @oma_genes = sort keys %{ $data_ref->{'uniprot'}->{$oma_uniprot}->{'gene'} };
            }
            else {
                warn "Cannot map OMA UniProt $oma_uniprot to gene in: $input\n";
            }
        }
        @oma_genes = sort @oma_genes;
        @oma_genes = uniq(@oma_genes);

        # Now, map from the *genes* to the OMA families.
        foreach my $oma_gene (@oma_genes) {
            $data_ref->{'gene'}->{$oma_gene}->{'oma_family'}->{$oma_family} = 1;
        }
    }
    else {
	die "In oma_groups file $oma_groups, cannot parse: $input\n";
    }
}
close $OMA_GROUPS;

my @oma_genes = sort keys %{ $data_ref->{'gene'} };
foreach my $oma_gene (@oma_genes) {
    my @oma_families = sort keys %{ $data_ref->{'gene'}->{$oma_gene}->{'oma_family'} };
    foreach my $oma_family (@oma_families) {
        print "$oma_gene\t$oma_family\n";
    }
}

