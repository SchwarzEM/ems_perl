#!/usr/bin/env perl

# sort_expr_levels_26feb2016.pl -- Erich Schwarz <ems394@cornell.edu>, 2/26/2016.
# Purpose: given gene annot file with well-headered TPM and minTPM, and known negative gene, generate lists and counts of genes expressed above background.

use strict;
use warnings;
use autodie;

use List::Util qw(max min);
use Scalar::Util qw(looks_like_number);

my $data_ref;

my $header1 = "Gene\tTPM\tminTPM\tGFP_tpm";
my $header2 = "Replicate\tExpr_tair_gene_count\tMax_TPM\tMin_TPM\tMax_minTPM\tMin_minTPM\tGFP_TPM";
my $header3 = "Genotype\tExpr_gene_count";

my $tair_gene_count = q{};

while (my $input = <>) {
    chomp $input;
    # from the header line, identify which columns contain which data types.
    if ( $input =~ /\A (\S+) \t (.+\S) \z/xms ) { 
        my $gene_id    = $1;
        my $field_text = $2;
        my $i          = 0;
        my @fields = split /\t/, $field_text;
        foreach my $field (@fields) {
            $i++;
            if ( ( $gene_id eq 'Gene' ) and ( $field =~ /\A (\S+_rep\d+) _ (TPM|minTPM) \z/xms ) ) {
                my $replicate = $1;
                my $data_type = $2;
                $data_ref->{'field'}->{$i}->{'replicate'} = $replicate;
                $data_ref->{'field'}->{$i}->{'data_type'} = $data_type;
            }
            elsif ( $gene_id ne 'Gene' ) { 
                if (     ( exists $data_ref->{'field'}->{$i}->{'replicate'} ) 
                     and ( exists $data_ref->{'field'}->{$i}->{'data_type'} ) ) { 

                    my $replicate = $data_ref->{'field'}->{$i}->{'replicate'};
                    my $data_type = $data_ref->{'field'}->{$i}->{'data_type'};

                    # Require values to be numerical.
                    if (looks_like_number($field) ) {
                        # Order the data in a way that makes it easy to mine efficiently later:
                        $data_ref->{'replicate'}->{$replicate}->{'gene'}->{$gene_id}->{'data_type'}->{$data_type} = $field;
                        
                    }
                }
            }
        }
    } 
    else { 
        die "Cannot parse input line: $input\n";
    }
}

my $replicate_summary = "sepal_replicate_expr_data.txt";
$replicate_summary    = safename($replicate_summary);   
open my $REPLICATE, '>', $replicate_summary;

my @replicates = sort keys %{ $data_ref->{'replicate'} };

foreach my $replicate (@replicates) {
    my $genotype = $replicate;
    if ( $genotype =~ /\A (\S+) _rep\d \z/xms ) { 
        $genotype = $1;
    }
    else {
        die "Cannot parse genotype from $genotype\n";
    }

    my @tair_genes      = grep { /AT(?:1|2|3|4|5|C|M)G\d+/ } keys %{ $data_ref->{'replicate'}->{$replicate}->{'gene'} };
    my $gfp_tpm         = $data_ref->{'replicate'}->{$replicate}->{'gene'}->{'GFP_coding_seq'}->{'data_type'}->{'TPM'};
    my @expr_tair_genes = ();
    $tair_gene_count   = @tair_genes;

    foreach my $tair_gene (@tair_genes) {
        my $min_tpm = $data_ref->{'replicate'}->{$replicate}->{'gene'}->{$tair_gene}->{'data_type'}->{'minTPM'};
        if ( $min_tpm > $gfp_tpm ) { 
            push @expr_tair_genes, $tair_gene;
        }
    }

    my $rep_gene_tpms = "$replicate.gene_tpms.tsv.txt";
    $rep_gene_tpms    = safename($rep_gene_tpms); 
    open my $TPM, '>', $rep_gene_tpms;

    @expr_tair_genes = sort @expr_tair_genes;
    my @expr_tair_gene_tpms     = ();
    my @expr_tair_gene_min_tpms = ();
    foreach my $expr_tair_gene (@expr_tair_genes) {
        my $expr_tair_gene_tpm     = $data_ref->{'replicate'}->{$replicate}->{'gene'}->{$expr_tair_gene}->{'data_type'}->{'TPM'};
        my $expr_tair_gene_min_tpm = $data_ref->{'replicate'}->{$replicate}->{'gene'}->{$expr_tair_gene}->{'data_type'}->{'minTPM'};

        print $TPM "$header1\n" if $header1;
        $header1 = q{};
        print $TPM "$expr_tair_gene\t$expr_tair_gene_tpm\t$expr_tair_gene_min_tpm\t$gfp_tpm\n";

        push @expr_tair_gene_tpms, $expr_tair_gene_tpm;
        push @expr_tair_gene_min_tpms, $expr_tair_gene_min_tpm;
        $data_ref->{'genotype'}->{$genotype}->{'expr_gene'}->{$expr_tair_gene} = 1;
        $data_ref->{'all_expr_gene'}->{$expr_tair_gene} = 1;
    }

    close $TPM;

    my $expr_tair_gene_count   = @expr_tair_genes;
    $expr_tair_gene_count      = commify($expr_tair_gene_count);

    my $max_expr_tair_gene_tpm = max(@expr_tair_gene_tpms);
    $max_expr_tair_gene_tpm    = commify($max_expr_tair_gene_tpm);

    my $min_expr_tair_gene_tpm = min(@expr_tair_gene_tpms);
    $min_expr_tair_gene_tpm    = commify($min_expr_tair_gene_tpm);

    my $max_expr_tair_gene_min_tpm = max(@expr_tair_gene_min_tpms);
    $max_expr_tair_gene_min_tpm    = commify($max_expr_tair_gene_min_tpm);

    my $min_expr_tair_gene_min_tpm = min(@expr_tair_gene_min_tpms);
    $min_expr_tair_gene_min_tpm    = commify($min_expr_tair_gene_min_tpm);

    print $REPLICATE "$header2\n" if $header2;
    $header2 = q{};
    print $REPLICATE "$replicate\t$expr_tair_gene_count\t",
                     "$max_expr_tair_gene_tpm\t$min_expr_tair_gene_tpm\t",
                     "$max_expr_tair_gene_min_tpm\t$min_expr_tair_gene_min_tpm\t",
                     "$gfp_tpm\n",
                     ;
}

close $REPLICATE;

my $overall_summary = "sepal_overall_expr_data.txt";
$overall_summary    = safename($overall_summary);
open my $OVERALL, '>', $overall_summary;

my @genotypes = sort keys %{ $data_ref->{'genotype'} };
foreach my $genotype (@genotypes) {
    my @expr_genes      = sort keys %{ $data_ref->{'genotype'}->{$genotype}->{'expr_gene'} };
    my $expr_gene_count = @expr_genes;
    $expr_gene_count    = commify($expr_gene_count);

    print $OVERALL "\n$header3\n" if $header3;
    $header3 = q{};
    print $OVERALL "$genotype\t$expr_gene_count\n";
}

my @all_expr_genes      = sort keys %{ $data_ref->{'all_expr_gene'} };
my $all_expr_gene_count = @all_expr_genes;
$all_expr_gene_count    = commify($all_expr_gene_count);
$tair_gene_count        = commify($tair_gene_count);

print $OVERALL "\nTotal expr. gene count: $all_expr_gene_count\n";
print $OVERALL "Total gene count: $tair_gene_count\n\n";

close $OVERALL;


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

