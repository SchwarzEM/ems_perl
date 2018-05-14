#!/usr/bin/env perl

# summarize_rsem_exp_01may2016.pl -- Erich Schwarz <ems394@cornell.edu>, 5/1/2016.
# Purpose: given gene annot file with well-headered TPM and minTPM, generate lists and counts of genes expressed above background.

use strict;
use warnings;
use autodie;

use List::Util qw(max min);
use Scalar::Util qw(looks_like_number);

my $data_ref;

my $header1 = "Gene\tCoding\tTPM\tminTPM";
my $header2 = "Gene\tCoding\tTPM\tminTPM";

my @header3_fields = qw(
    Bio_type
    Expr_conf_prot_gene_count
    Expr_conf_ncRNA_gene_count
    Expr_any_prot_gene_count
    Expr_any_ncRNA_gene_count
    Max_expr_conf_gene_TPM
    Max_expr_all_gene_TPM
    Min_expr_conf_gene_TPM
    Min_expr_all_gene_TPM
    Max_expr_conf_gene_minTPM
    Max_expr_any_gene_minTPM
    Min_expr_conf_gene_minTPM
    Min_expr_all_gene_minTPM
);
my $header3 = join "\t", @header3_fields;

my $gene_count = q{};

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
            if ( $gene_id eq 'Gene' ) { 
                if ( $field =~ /\A (\S+) _ (TPM|minTPM) \z/xms ) {
                    my $bio_type = $1;
                    my $data_type = $2;
                    $data_ref->{'field'}->{$i}->{'bio_type'}  = $bio_type;
                    $data_ref->{'field'}->{$i}->{'data_type'} = $data_type;
                }
                # Can generalize this test of $field, as needed:
                elsif ( $field eq 'Coding' ) {
                    $data_ref->{'field'}->{$i}->{'annot_type'} = $field;
                }
            }
            elsif ( $gene_id ne 'Gene' ) { 
                if (     ( exists $data_ref->{'field'}->{$i}->{'bio_type'} ) 
                     and ( exists $data_ref->{'field'}->{$i}->{'data_type'} ) ) { 

                    my $bio_type  = $data_ref->{'field'}->{$i}->{'bio_type'};
                    my $data_type = $data_ref->{'field'}->{$i}->{'data_type'};

                    # Require values to be numerical.
                    if (looks_like_number($field) ) {
                        # Order the data in a way that makes it easy to mine efficiently later:
                        $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$gene_id}->{'data_type'}->{$data_type} = $field;
                    }
                }
                elsif ( $data_ref->{'field'}->{$i}->{'annot_type'} ) {
                    my $annot_type = $data_ref->{'field'}->{$i}->{'annot_type'};
                     $data_ref->{'gene'}->{$gene_id}->{'annot_type'}->{$annot_type} = $field;
                }
            }
        }
    } 
    else { 
        die "Cannot parse input line: $input\n";
    }
}

my $bio_type_summary = 'bio_type_expr_data.txt';
$bio_type_summary    = safename($bio_type_summary);   
open my $BIO_TYPE, '>', $bio_type_summary;

my @bio_types = sort keys %{ $data_ref->{'bio_type'} };

foreach my $bio_type (@bio_types) {
    my @genes                 = sort keys %{ $data_ref->{'bio_type'}->{$bio_type}->{'gene'} };
    my @expr_prot_genes_conf  = ();
    my @expr_ncrna_genes_conf = ();
    my @expr_prot_genes_any   = ();
    my @expr_ncrna_genes_any  = ();
    $gene_count               = @genes;

    foreach my $gene (@genes) {
        my $tpm     = $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$gene}->{'data_type'}->{'TPM'};
        my $min_tpm = $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$gene}->{'data_type'}->{'minTPM'};
        my $coding  = $data_ref->{'gene'}->{$gene}->{'annot_type'}->{'Coding'};
        if ( $tpm     >= 0.1 ) {
            if ( $coding =~ /protein/xms ) { 
                push @expr_prot_genes_any, $gene;
            }
            elsif ( $coding eq 'ncRNA' ) {
                push @expr_ncrna_genes_any, $gene;
            }
            else {
                die "Cannot decipher coding status of gene \"$gene\": $coding\n";
            }
        }
        if ( $min_tpm >= 0.1 ) { 
            if ( $coding =~ /protein/xms ) {
                push @expr_prot_genes_conf, $gene;
            }
            elsif ( $coding eq 'ncRNA' ) {
                push @expr_ncrna_genes_conf, $gene;
            }       
            else {
                die "Cannot decipher coding status of gene \"$gene\": $coding\n";
            }
        }
    }

    my $rep_gene_conf_tpms = "$bio_type.conf.gene_tpms.tsv.txt";
    my $rep_gene_any_tpms  = "$bio_type.any.gene_tpms.tsv.txt";
    $rep_gene_conf_tpms    = safename($rep_gene_conf_tpms); 
    $rep_gene_any_tpms     = safename($rep_gene_any_tpms);

    @expr_prot_genes_conf  = sort @expr_prot_genes_conf;
    @expr_ncrna_genes_conf = sort @expr_ncrna_genes_conf;
    my @expr_genes_conf    = (@expr_prot_genes_conf,@expr_ncrna_genes_conf);
    @expr_genes_conf       = sort @expr_genes_conf;

    @expr_prot_genes_any  = sort @expr_prot_genes_any;
    @expr_ncrna_genes_any = sort @expr_ncrna_genes_any;
    my @expr_genes_any    = (@expr_prot_genes_any,@expr_ncrna_genes_any);
    @expr_genes_any       = sort @expr_genes_any;

    my @expr_gene_conf_tpms     = ();
    my @expr_gene_conf_min_tpms = ();

    my @expr_gene_any_tpms     = ();
    my @expr_gene_any_min_tpms = ();

    open my $TPM_CONF, '>', $rep_gene_conf_tpms;
    foreach my $expr_gene_conf (@expr_genes_conf) {
        my $expr_gene_tpm     = $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$expr_gene_conf}->{'data_type'}->{'TPM'};
        my $expr_gene_min_tpm = $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$expr_gene_conf}->{'data_type'}->{'minTPM'};
        my $coding            = $data_ref->{'gene'}->{$expr_gene_conf}->{'annot_type'}->{'Coding'};

        print $TPM_CONF "$header1\n" if $header1;
        $header1 = q{};
        print $TPM_CONF "$expr_gene_conf\t$coding\t$expr_gene_tpm\t$expr_gene_min_tpm\n";

        push @expr_gene_conf_tpms,     $expr_gene_tpm;
        push @expr_gene_conf_min_tpms, $expr_gene_min_tpm;
        $data_ref->{'bio_type'}->{$bio_type}->{'expr_conf_gene'}->{$expr_gene_conf} = 1;
    }
    close $TPM_CONF;

    open my $TPM_ANY,  '>', $rep_gene_any_tpms;
    foreach my $expr_gene_any (@expr_genes_any) {
        my $expr_gene_tpm     = $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$expr_gene_any}->{'data_type'}->{'TPM'};
        my $expr_gene_min_tpm = $data_ref->{'bio_type'}->{$bio_type}->{'gene'}->{$expr_gene_any}->{'data_type'}->{'minTPM'};
        my $coding            = $data_ref->{'gene'}->{$expr_gene_any}->{'annot_type'}->{'Coding'};

        print $TPM_ANY "$header2\n" if $header2;
        $header2 = q{};
        print $TPM_ANY "$expr_gene_any\t$coding\t$expr_gene_tpm\t$expr_gene_min_tpm\n";

        push @expr_gene_any_tpms,     $expr_gene_tpm;
        push @expr_gene_any_min_tpms, $expr_gene_min_tpm;
        $data_ref->{'bio_type'}->{$bio_type}->{'expr_any_gene'}->{$expr_gene_any} = 1;
    }
    close $TPM_ANY;

    my $expr_conf_prot_gene_count  = @expr_prot_genes_conf;
    my $expr_conf_ncrna_gene_count = @expr_ncrna_genes_conf;
    $expr_conf_prot_gene_count     = commify($expr_conf_prot_gene_count);
    $expr_conf_ncrna_gene_count    = commify($expr_conf_ncrna_gene_count);

    my $expr_any_prot_gene_count  = @expr_prot_genes_any;
    my $expr_any_ncrna_gene_count = @expr_ncrna_genes_any;
    $expr_any_prot_gene_count     = commify($expr_any_prot_gene_count);
    $expr_any_ncrna_gene_count    = commify($expr_any_ncrna_gene_count);

    my $max_expr_conf_gene_tpm = max(@expr_gene_conf_tpms);
    $max_expr_conf_gene_tpm    = commify($max_expr_conf_gene_tpm);

    my $max_expr_all_gene_tpm = max(@expr_gene_any_tpms);
    $max_expr_all_gene_tpm    = commify($max_expr_all_gene_tpm);

    my $min_expr_conf_gene_tpm = min(@expr_gene_conf_tpms);
    $min_expr_conf_gene_tpm    = commify($min_expr_conf_gene_tpm);

    my $min_expr_all_gene_tpm = min(@expr_gene_any_tpms);
    $min_expr_all_gene_tpm    = commify($min_expr_all_gene_tpm);

    my $max_expr_conf_gene_min_tpm = max(@expr_gene_conf_min_tpms);
    $max_expr_conf_gene_min_tpm    = commify($max_expr_conf_gene_min_tpm);

    my $max_expr_any_gene_min_tpm = max(@expr_gene_any_min_tpms);
    $max_expr_any_gene_min_tpm    = commify($max_expr_any_gene_min_tpm);

    my $min_expr_conf_gene_min_tpm = min(@expr_gene_conf_min_tpms);
    $min_expr_conf_gene_min_tpm    = commify($min_expr_conf_gene_min_tpm);

    my $min_expr_all_gene_min_tpm = min(@expr_gene_any_min_tpms);
    $min_expr_all_gene_min_tpm    = commify($min_expr_all_gene_min_tpm);

    print $BIO_TYPE "$header3\n" if $header3;
    $header3 = q{};

    my $output = $bio_type
                 . "\t"
                 . $expr_conf_prot_gene_count
                 . "\t"
                 . $expr_conf_ncrna_gene_count
                 . "\t"
                 . $expr_any_prot_gene_count
                 . "\t"
                 . $expr_any_ncrna_gene_count
                 . "\t"
                 . $max_expr_conf_gene_tpm
                 . "\t"
                 . $max_expr_all_gene_tpm
                 . "\t"
                 . $min_expr_conf_gene_tpm
                 . "\t"
                 . $min_expr_all_gene_tpm
                 . "\t"
                 . $max_expr_conf_gene_min_tpm
                 . "\t"
                 . $max_expr_any_gene_min_tpm
                 . "\t"
                 . $min_expr_conf_gene_min_tpm
                 . "\t"
                 . $min_expr_all_gene_min_tpm
                 ;

    print $BIO_TYPE "$output\n";
}
close $BIO_TYPE;

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

