#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Statistics::Descriptive;

my $genes        = 0;

my @prot_sizes   = ();

my $prots_100aa  = 0;
my $prots_99aa   = 0;

my $phobius_sig_count    = 0;
my $phobius_tm_count     = 0;
my $phobius_sig_tm_count = 0;

my $pfam_count   = 0;

my @obs_tpms    = ();
my $max_obs_tpm = 0;
my @tpms        = ();
my $expr_genes  = 0;

my @e_values     = ();
my @geo_e_values = ();

# Headers of data:
#
# Gene	Prot_size	Max_prot_size	Phobius	NCoils	Psegs	PFAM	OFind_4spp	OFind_4spp_Summary	OFind_5spp	
# OFind_5spp_Summary	OFind_6spp	OFind_6spp_Summary	OFind_8spp	OFind_8spp_Summary	Male.v.Fem.logFC	
# Male.vs.Fem.FDR	Male.v.modENC.logFC	Male.vs.modENC.FDR	Fem.v.modENC.logFC	Fem.vs.modENC.FDR	Male_1_TPM	
# Male_1_reads	Male_2_TPM	Male_2_reads	Male_3_TPM	Male_3_reads	Female_1_TPM	Female_1_reads	Female_2_TPM	Female_2_reads	
# Female_3_TPM	Female_3_reads	modENCODE_TPM	modENCODE_reads	E_value_BlastP_5spp

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t   # $1 = Gene
                       [^\t]* \t
                       (\S+) \t      # $2 = Max_prot_size
                       ([^\t]*) \t   # $3 = Phobius
                       (?: [^\t]* \t){2}
                       ([^\t]*) \t   # $4 = PFAM
                       (?: [^\t]* \t){8}
                       ([^\t]*) \t   # $5 = Male.v.Fem.logFC
                       ([^\t]*) \t   # $6 = Male.vs.Fem.FDR
                       (?: [^\t]* \t){4}
                       ([^\t]*) \t   # $7 = Male_1_TPM
                       [^\t]* \t
                       ([^\t]*) \t   # $8 = Male_2_TPM
                       [^\t]* \t
                       ([^\t]*) \t   # $9 = Male_3_TPM
                       [^\t]* \t
                       ([^\t]*) \t   # $10 = Female_1_TPM  
                       [^\t]* \t
                       ([^\t]*) \t   # $11 = Female_2_TPM
                       [^\t]* \t
                       ([^\t]*) \t   # $12 = Female_3_TPM
                       [^\t]* \t
                       ([^\t]*) \t   # $13 = modENCODE_TPM
                       [^\t]* \t
                       ([^\t]*)      # $14 = E_value_BlastP_5spp
    \z/xms ) {
        my $gene = $1;
        my $max_prot_size = $2;
        my $phobius = $3;
        my $pfam = $4;
        my $male_v_fem_log_fc = $5;
        my $male_vs_fem_fdr = $6;
        my $male_1_tpm = $7;
        my $male_2_tpm = $8;
        my $male_3_tpm = $9;
        my $female_1_tpm = $10;
        my $female_2_tpm = $11;
        my $female_3_tpm = $12;
        my $modencode_tpm = $13;
        my $e_value_blastp_5spp = $14;

        @obs_tpms    = ();
        $max_obs_tpm = 0;

        if ( $gene ne 'Gene' ) { 
            $genes++;  # basic count of data lines

            push @prot_sizes, $max_prot_size;

            if ( $max_prot_size >= 100 ) {
                $prots_100aa++;
            }
            else {
                $prots_99aa++;
            }

            if ( ( $phobius =~ /SigP/xms ) and ( $phobius =~ /TM/xms ) ) {
                $phobius_sig_tm_count++;
            }
            elsif ( $phobius =~ /SigP/xms ) {
                $phobius_sig_count++;
            }
            elsif ( $phobius =~ /TM/xms ) {
                $phobius_tm_count++;
            }

            if ( $pfam =~ / PF\d+ /xms ) { 
                 $pfam_count++;
            }

            @obs_tpms = ($male_1_tpm, $male_2_tpm, $male_3_tpm, $female_1_tpm, $female_2_tpm, $female_3_tpm, $modencode_tpm);
            $max_obs_tpm = max(@obs_tpms);
            if ( $max_obs_tpm >= 0.1 ) {
                $expr_genes++;
            }
            push @tpms, $max_obs_tpm;

            # Only count E-values that actually exist.
            if ( $e_value_blastp_5spp =~ /\S/xms ) { 
                push @e_values, $e_value_blastp_5spp;
                # Can't do geometric means with '0', so set set arbitrary minimum of 1e-200.
                if ( $e_value_blastp_5spp == 0 ) {
                    $e_value_blastp_5spp = 1e-200;
                }
                push @geo_e_values, $e_value_blastp_5spp;
            }
        }
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}

         
my $stat_prot_sizes = Statistics::Descriptive::Full->new();
$stat_prot_sizes->add_data(@prot_sizes);

my $mean_prot_sizes   = $stat_prot_sizes->mean();
my $median_prot_sizes = $stat_prot_sizes->median();
my $max_prot_sizes    = $stat_prot_sizes->max();
my $min_prot_sizes    = $stat_prot_sizes->min();

my $stat_tpms = Statistics::Descriptive::Full->new();
$stat_tpms->add_data(@tpms);

my $mean_tpms = $stat_tpms->mean();
my $median_tpms = $stat_tpms->median();
my $max_tpms = $stat_tpms->max();
my $min_tpms = $stat_tpms->min();

my $stat_e_values = Statistics::Descriptive::Full->new();
$stat_e_values->add_data(@e_values);

my $stat_geo_e_values = Statistics::Descriptive::Full->new();
$stat_geo_e_values->add_data(@geo_e_values);

my $count_e_values    = $stat_e_values->count();
my $geo_mean_e_values = $stat_geo_e_values->geometric_mean();
my $median_e_values   = $stat_e_values->median();
my $max_e_values      = $stat_e_values->max();
my $min_e_values      = $stat_e_values->min();

print "\n";

my $genes_commmaed = commify($genes);
print "Number of genes:          $genes_commmaed\n";

print "\n";

my $ratio = ($prots_100aa/$genes);
$ratio = sprintf("%.3f", $ratio);
$prots_100aa = commify($prots_100aa);
print "Prots. >= 100 aa:         $prots_100aa ($ratio)\n";

$ratio = ($prots_99aa/$genes);
$ratio = sprintf("%.3f", $ratio);
$prots_99aa = commify($prots_99aa);
print "Prots. <= 99 aa:          $prots_99aa ($ratio)\n";

$mean_prot_sizes = commify($mean_prot_sizes);
$mean_prot_sizes = sprintf("%.1f", $mean_prot_sizes);
print "Mean protein size:        $mean_prot_sizes\n";

$median_prot_sizes = commify($median_prot_sizes);
print "Median protein size:      $median_prot_sizes\n";

$max_prot_sizes = commify($max_prot_sizes);
print "Maximum protein size:     $max_prot_sizes\n";

$min_prot_sizes = commify($min_prot_sizes);
print "Minimum protein size:     $min_prot_sizes\n";

print "\n";

$ratio = ($phobius_sig_count/$genes);
$ratio = sprintf("%.3f", $ratio);
$phobius_sig_count = commify($phobius_sig_count);
print "No. of SigP proteins:     $phobius_sig_count ($ratio)\n";

$ratio = ($phobius_tm_count/$genes);
$ratio = sprintf("%.3f", $ratio);
$phobius_tm_count = commify($phobius_tm_count);
print "No. of TM proteins:       $phobius_tm_count ($ratio)\n";

$ratio = ($phobius_sig_tm_count/$genes);
$ratio = sprintf("%.3f", $ratio);
$phobius_sig_tm_count = commify($phobius_sig_tm_count);
print "No. of SigP+TM proteins:  $phobius_sig_tm_count ($ratio)\n";

print "\n";

$ratio = ($pfam_count/$genes);
$ratio = sprintf("%.3f", $ratio);
$pfam_count = commify($pfam_count);
print "Genes enc. PFAM domain:   $pfam_count ($ratio)\n";

print "\n";

$ratio = ($expr_genes/$genes);
$ratio = sprintf("%.3f", $ratio);
$expr_genes = commify($expr_genes);
print "Expr. >=0.1 TPM:  $expr_genes ($ratio)\n";

print "\n";

$mean_tpms = sprintf("%.2f", $mean_tpms);
$mean_tpms = commify($mean_tpms);
print "Mean TPM:    $mean_tpms\n";

$median_tpms = sprintf("%.2f", $median_tpms);
$median_tpms = commify($median_tpms);
print "Median TPM:  $median_tpms\n";

$max_tpms = sprintf("%.2f", $max_tpms);
$max_tpms = commify($max_tpms);
print "Max. TPM:    $max_tpms\n";

$min_tpms = sprintf("%.2f", $min_tpms);
$min_tpms = commify($min_tpms);
print "Min. TPM:    $min_tpms\n";

print "\n";

$ratio = ($count_e_values/$genes);
$ratio = sprintf("%.3f", $ratio);
$count_e_values = commify($count_e_values);
print "Genes w/ BlastP E-values:    $count_e_values ($ratio)\n";

print "\n";

print "Mean (geo.) BlastP E-value:  $geo_mean_e_values\n";
print "Median BlastP E-value:       $median_e_values\n";
print "Maximum BlastP E-value:      $max_e_values\n";
print "Minimum BlastP E-value:      $min_e_values\n";
print "\n";


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

