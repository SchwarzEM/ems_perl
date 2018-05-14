#!/usr/bin/env perl

# map_gene2ref_density.pl -- Erich Schwarz <ems394@cornell.edu>, 8/21/2013:
# Purpose: given a UniProt to gene table, and given a UniProtKB file mapping UniProts to primary refs. and electronically mapped refs., get reference density/gene.

# Note that this program weights references in parallel: primary refs. are kept separate from electronically mapped ones, 
#    because a primary ref in UniProt by definition has been manually annotated, and can be assumed to carry more specific information
#    about a primary gene even if it has also been electronically mapped to N other genes.

# Also note that this uses a processed version of the UniProt-to-PubMed ID table, not UniProtKB's raw file.
# To get the former from the latter, try:  cut -f 1,5,17 HUMAN_9606_idmapping_selected.tab > processed.txt.

use strict;
use warnings;

use List::Util qw(sum);
use List::MoreUtils qw(uniq);

my $data_ref;

my $uniprot2gene = $ARGV[0];
my $uniprot2ref  = $ARGV[1];

open my $PROT2GENE, '<', $uniprot2gene or die "Can't open UniProt to gene table $uniprot2gene: $!";
while (my $input = <$PROT2GENE>) { 
    chomp $input;
    # N.B. don't require that there be *only* two fields (so an input file with descriptions in the third field is OK).
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $uniprot = $1;
        my $gene    = $2;

        # Note that UniProt allows one protein to come from two or more genes; ergo, mapping is not many-to-unique.
        $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;

        # Recording this ensures that we list *all* genes, not just ones that have references.
        $data_ref->{'gene'}->{$gene}->{'seen'} = 1;
    }
    else { 
        die "From UniProt to gene table $uniprot2gene, can't parse input line: $input\n";
    }
}
close $PROT2GENE or die "Can't close filehandle to UniProt to gene table $uniprot2gene: $!";

open my $PROT2REF, '<', $uniprot2ref or die "Can't open UniProt to PubMed ID ref. table $uniprot2ref: $!";
while (my $input = <$PROT2REF>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]*) \t ([^\t]*) /xms ) { 
        my $uniprot      = $1;
        my $prim_ref_txt = $2;
        my $sec_ref_txt  = $3;
        my @prim_refs = split /;\s+/, $prim_ref_txt;
        my @sec_refs  = split /;\s+/, $sec_ref_txt;
        if ( exists $data_ref->{'uniprot'}->{$uniprot}->{'gene'} ) { 
            my @genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene'} };
            foreach my $gene (@genes) { 
                foreach my $prim_ref (@prim_refs) { 
                    $data_ref->{'gene'}->{$gene}->{'prim_ref'}->{$prim_ref} = 1;
                    # Also map the other way, so that we can get citation densities.
                    $data_ref->{'prim_ref'}->{$prim_ref}->{'gene'}->{$gene} = 1;
                }
                foreach my $sec_ref (@sec_refs) {
                    $data_ref->{'gene'}->{$gene}->{'sec_ref'}->{$sec_ref} = 1;
                    $data_ref->{'sec_ref'}->{$sec_ref}->{'gene'}->{$gene} = 1;
                }
            }
        }
    }
    else {
        die "From UniProt to PubMed ID ref. table $uniprot2ref, can't parse input line: $input\n";
    }     
}
close $PROT2REF or die "Can't close filehandle to UniProt to PubMed ID ref. table $uniprot2ref: $!";

# Compute citation densities.
my @prim_refs1 = sort keys %{ $data_ref->{'prim_ref'} };
foreach my $prim_ref1 (@prim_refs1) {
    my @genes            = sort keys %{ $data_ref->{'prim_ref'}->{$prim_ref1}->{'gene'} };
    my $gene_count       = @genes;

    my $citation_density = (1 / $gene_count);
    # To lower floating-point precision to slightly less absurd levels:
    $citation_density    = sprintf "%.4f", $citation_density;

    $data_ref->{'prim_ref'}->{$prim_ref1}->{'value'} = $citation_density;
}
my @sec_refs1 = sort keys %{ $data_ref->{'sec_ref'} };
foreach my $sec_ref1 (@sec_refs1) {
    my @genes            = sort keys %{ $data_ref->{'sec_ref'}->{$sec_ref1}->{'gene'} };
    my $gene_count       = @genes;

    my $citation_density = (1 / $gene_count);
    # To lower floating-point precision to slightly less absurd levels:
    $citation_density    = sprintf "%.4f", $citation_density;

    $data_ref->{'sec_ref'}->{$sec_ref1}->{'value'} = $citation_density;
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) { 
    my @total_refs      = ();

    my $total_ref_score = 0;
    my $prim_ref_score  = 0;
    if ( exists $data_ref->{'gene'}->{$gene}->{'prim_ref'} ) { 
        my @prim_refs   = sort keys %{ $data_ref->{'gene'}->{$gene}->{'prim_ref'} };
        my @prim_values = map { $data_ref->{'prim_ref'}->{$_}->{'value'} } @prim_refs;
        $prim_ref_score = sum @prim_values;
        push @total_refs, @prim_refs;
    }
    my $sec_ref_score = 0;
    if ( exists $data_ref->{'gene'}->{$gene}->{'sec_ref'} ) {
        my @sec_refs   = sort keys %{ $data_ref->{'gene'}->{$gene}->{'sec_ref'} };  
        my @sec_values = map { $data_ref->{'sec_ref'}->{$_}->{'value'} } @sec_refs;
        $sec_ref_score = sum @sec_values;
        push @total_refs, @sec_refs;
    }
    $total_ref_score = ($prim_ref_score + $sec_ref_score);

    @total_refs     = sort @total_refs;
    @total_refs     = uniq @total_refs;
    my $ref_summary = join '; ', @total_refs;

    print "$gene\t$total_ref_score\t$ref_summary\n";
}


