#!/usr/bin/env perl

# make_GO_term_matrix_18mar2013a.pl -- Erich Schwarz <ems394@cornell.edu>, 3/19/2013.
# Purpose: given many refined FUNC outputs, generate a matrix in which *GO terms* are the genes, up- or down-regulations are the stages, and the values are either 0/1 (for absence/presence of an enriched GO term) or 0/N (for absence/presence, where N = -log10[p-value]).

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(zip);

my $binary;
my $log_scale;
my @ontologies    = ();
my %ok_ontologies = ( 'molecular_function' => 1,
                      'cellular_component' => 1,
                      'biological_process' => 1, );
my $debug;
my $help;

GetOptions ( 'binary'           => \$binary,
             'log_scale'        => \$log_scale,
             'ontologies=s{1,}' => \@ontologies,
             'debug'            => \$debug,
             'help'             => \$help,  );

if ( $help or ( (! $binary ) and (! $log_scale ) ) or ( $binary and $log_scale ) ) {
    die "Format: make_GO_term_matrix_18mar2013a.pl\n",
        "            --log_scale|-l [or] --binary [to get N/0 vs. 1/0 matrices; required, and mutually exclusive options]\n",
        "            --ontologies|-o [specify 1+ of: \"molecular_function cellular_component biological_process\"; use to get subset matrices; default is all three]\n",
        "            --debug|-d [print minimum observed non-zero p-values to STDERR]\n",
        "            --help|-h [print this message]\n",
        ;
}

if (! @ontologies) {
    @ontologies = qw( molecular_function cellular_component biological_process );
}

foreach my $ontology (@ontologies) {
    if (! exists $ok_ontologies{$ontology} ) { 
        die "Unacceptable ontology: $ontology\n";
    }
}

my @data_sets = qw( acey_dirs cel_b2go_dirs cel_wbGO_dirs hco_Pasi_dirs hco_aug_dirs );

my @condition_up_list   = ();
my @condition_down_list = ();
my @full_condition_list = ();

my %data_synonym = ( 'acey_dirs'     => "Acey_", 
                     'cel_b2go_dirs' => "Cel_",
                     'cel_wbGO_dirs' => "Cel_WS230_",
                     'hco_Pasi_dirs' => "Hco_",
                     'hco_aug_dirs'  => "Hco_aug_", );

my $data_ref;

@{ $data_ref->{'acey_dirs'} } = qw( 
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.24.PI.vs.L3i.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.24HCM.vs.L3i.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.24.PI.vs.24HCM.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.5.D.vs.24.PI.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.12.D.vs.5.D.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.17.D.vs.12.D.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.19.D.vs.17.D.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.Alb.4hr.D18.vs.noAlb.cont.D18.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.18D.Cry5B.4hr.vs.18D.HEPES.4hr.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.18D.Cry5B.24hr.vs.18D.HEPES.24hr.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Acey_RSEM_03dec2012.log10.18D.SB.plusCry5B.vs.18D.SB.plusHEPES.func_wilcoxon_outdir 
);

@{ $data_ref->{'cel_b2go_dirs'} } = qw ( 
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.late_embryo.vs.early_embryo.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.b2go.log10.L1.vs.late_embryo.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.embryo_30min.vs.embryo_0min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.embryo_60min.vs.embryo_30min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.embryo_90min.vs.embryo_60min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.embryo_120min.vs.embryo_90min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.b2go.log10.L1.vs.embryo_120min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.b2go.log10.L2.vs.L1.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.L3.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.herm_L4.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.b2go.log10.male_L4.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.b2go.log10.male_L4.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.YA.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17mar2013.b2go.log10.dauer_entry.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.dauer_entry.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.dauer_entry.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.dauer_exit.vs.dauer_entry.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.dauer_exit.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.dauer_exit.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.YA.vs.dauer_exit.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.Alb.4hr.L4.vs.noAlb.cont.L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.Alb.4hr.L4.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.noAlb.cont.L4.vs.herm_L4.func_wilcoxon_outdir 
);

@{ $data_ref->{'cel_wbGO_dirs'} } = qw ( 
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.late_embryo.vs.early_embryo.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.wbGO.log10.L1.vs.late_embryo.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.embryo_30min.vs.embryo_0min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.embryo_60min.vs.embryo_30min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.embryo_90min.vs.embryo_60min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.embryo_120min.vs.embryo_90min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.wbGO.log10.L1.vs.embryo_120min.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.wbGO.log10.L2.vs.L1.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.L3.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.herm_L4.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.wbGO.log10.male_L4.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.wbGO.log10.male_L4.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.YA.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17mar2013.wbGO.log10.dauer_entry.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.dauer_entry.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.dauer_entry.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.dauer_exit.vs.dauer_entry.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.dauer_exit.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.dauer_exit.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.YA.vs.dauer_exit.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.Alb.4hr.L4.vs.noAlb.cont.L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.Alb.4hr.L4.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_05dec2012.log10.noAlb.cont.L4.vs.herm_L4.func_wilcoxon_outdir
);

@{ $data_ref->{'hco_Pasi_dirs'} } = qw ( 
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.L1.vs.Egg.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.L2.vs.L1.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.L3.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.L4.F.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.L4.M.vs.L3.func_wilcoxon_outdir  
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_18mar2013.log10.L4.M.vs.L4.F.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.Adult.F.vs.L4.F.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_11dec2012.log10.Adult.M.vs.L4.M.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.Pasi_RSEM_18mar2013.log10.Adult.M.vs.Adult.F.func_wilcoxon_outdir
);

@{ $data_ref->{'hco_aug_dirs'} } = qw ( 
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.L1.vs.Egg.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.L2.vs.L1.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.L3.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.L4.F.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.L4.M.vs.L3.func_wilcoxon_outdir  
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_18mar2013.log10.L4.M.vs.L4.F.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.Adult.F.vs.L4.F.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_11dec2012.log10.Adult.M.vs.L4.M.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Hco.aug_RSEM_18mar2013.log10.Adult.M.vs.Adult.F.func_wilcoxon_outdir
);

foreach my $data_set (@data_sets) {
    foreach my $data_dir ( @{ $data_ref->{$data_set} } ) {
        my $dir_basename = basename $data_dir;
        # Sample $dir_basename:
        # Acey_RSEM_03dec2012.log10.24.PI.vs.L3i.func_wilcoxon_outdir
        if ( $dir_basename =~ /\A \S+ _RSEM \S+ \.log10\. (\S+) \.vs\. (\S+) \.func_wilcoxon_outdir \z/xms ) { 
            my $numerator   = $1;
            my $denominator = $2;
            my $condition   = q{};

            my $upreg_gene_command = 'filter_func_refinement.pl'
                                     . " -i $data_dir"
                                     . '/refinement-*-0.001_0.001.txt -u | cut -f 1-3,8'
                                     ;

            my @upreg_term_text    = grep { $_ !~ /p_high_ranks_after_refinement/ } `$upreg_gene_command`;
            chomp @upreg_term_text;

            $condition = $data_synonym{$data_set} . "$denominator.to.$numerator" . '[+]';
            push @condition_up_list, $condition;

            foreach my $upreg_term (@upreg_term_text) { 
                if ( $upreg_term =~ /\A ([^\t]+) \t ([^\t]+) \t ([^\t]+) \t (\S+) \z/xms ) {
                    my $ontology = $1;
                    my $go_text  = $2;
                    my $go_id    = $3;
                    my $p_value  = $4;
                    my $go_term  = $go_id . "[$go_text]";
                    $go_term     =~ s/\s/_/g;

                    $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'condition'}->{$condition}->{'p-value'} = $p_value;

                    if ( exists $data_ref->{'min_obs_nonzero_p-value'} ) {
                        if ( ( $data_ref->{'min_obs_nonzero_p-value'} > $p_value ) and ( $p_value > 0 ) ) { 
                            $data_ref->{'min_obs_nonzero_p-value'} = $p_value;
                            if ($debug) {
                                warn "DEBUG: revised min_obs_nonzero_p-value downward: $data_ref->{'min_obs_nonzero_p-value'}\n";
                            }
                        }
                    }

                    # Note that if " $p_value > 0" is not made an explicit condition, a p-value of 0 can sneak in here, and then will 'exist'!
                    if ( (! exists $data_ref->{'min_obs_nonzero_p-value'} ) and ( $p_value > 0 ) ) { 
                        $data_ref->{'min_obs_nonzero_p-value'} = $p_value;
                        if ($debug) { 
                            warn "DEBUG: initialized min_obs_nonzero_p-value: $data_ref->{'min_obs_nonzero_p-value'}\n";
                        }
                    }
                }
                else {
                    die "Can't parse upregulated term $upreg_term\n";
                }
            }

            my $downreg_gene_command = 'filter_func_refinement.pl'
                                     . " -i $data_dir"
                                     . '/refinement-*-0.001_0.001.txt -d | cut -f 1-3,7'
                                     ;

            my @downreg_term_text  = grep { $_ !~ /p_low_ranks_after_refinement/ } `$downreg_gene_command`;
            chomp @downreg_term_text;

            $condition = $data_synonym{$data_set} . "$denominator.to.$numerator" . '[-]';
            push @condition_down_list, $condition;

            foreach my $downreg_term (@downreg_term_text) {
                if ( $downreg_term =~ /\A ([^\t]+) \t ([^\t]+) \t ([^\t]+) \t (\S+) \z/xms ) {
                    my $ontology = $1;
                    my $go_text  = $2;
                    my $go_id    = $3;
                    my $p_value  = $4;
                    my $go_term  = $go_id . "[$go_text]";
                    $go_term     =~ s/\s/_/g;

                    $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'condition'}->{$condition}->{'p-value'} = $p_value;

                    if ( exists $data_ref->{'min_obs_nonzero_p-value'} ) {
                        if ( ( $data_ref->{'min_obs_nonzero_p-value'} > $p_value ) and ( $p_value > 0 ) ) {
                            $data_ref->{'min_obs_nonzero_p-value'} = $p_value;
                            if ($debug) {
                                warn "DEBUG: revised min_obs_nonzero_p-value downward: $data_ref->{'min_obs_nonzero_p-value'}\n";
                            }
                        }
                    }
            
                    # Note that if " $p_value > 0" is not made an explicit condition, a p-value of 0 can sneak in here, and then will 'exist'!
                    if ( (! exists $data_ref->{'min_obs_nonzero_p-value'} ) and ( $p_value > 0 ) ) {
                        $data_ref->{'min_obs_nonzero_p-value'} = $p_value;
                        if ($debug) {
                            warn "DEBUG: initialized min_obs_nonzero_p-value: $data_ref->{'min_obs_nonzero_p-value'}\n";
                        }
                    }
                }
                else {
                    die "Can't parse downregulated term $downreg_term\n";
                }
            }
        }
        else {
            die "Can't parse specific label in $dir_basename\n";
        }
    }
}

@full_condition_list = zip(@condition_up_list,@condition_down_list);

my $header = join "\t", @full_condition_list;
$header    = "GO_term\t$header";
print "$header\n";

foreach my $ontology (@ontologies) { 
    my @go_terms = sort keys %{ $data_ref->{'ontology'}->{$ontology}->{'go_term'} };
    foreach my $go_term (@go_terms) { 
        my @text_fields = ();
        push @text_fields, $go_term;
        foreach my $condition (@full_condition_list) {
            my $p_value           = q{};
            my $neg_log10_p_value = q{};

            # If a p-value was recorded, set it to min_obs_nonzero_p-value, then get its -log10.
            if ( exists $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'condition'}->{$condition}->{'p-value'} ) { 
                $p_value = $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'condition'}->{$condition}->{'p-value'};
                if ( ( $p_value < 0 ) or ( $p_value > 1 ) ) {
                    die "Error -- impossible 'p-value' of $p_value\n";
                }
                else {
                    if ( $p_value < $data_ref->{'min_obs_nonzero_p-value'} ) {
                        $p_value = $data_ref->{'min_obs_nonzero_p-value'} ;
                    }
                    $neg_log10_p_value = ( -1 * ( log($p_value)/log(10) ) );
                }
            }

            # If a p-value was *not* recorded, assume that its p-value is 1, and that its -log10 is 0.
            else { 
                $neg_log10_p_value = 0;
            }

            # If we're doing just 'binary', make everything very simple:
            if ( $binary and ( $neg_log10_p_value > 0 ) ) { 
                $neg_log10_p_value = 1;
            }

            push @text_fields, $neg_log10_p_value;
        }
        my $output_text = join "\t", @text_fields;
        print "$output_text\n";
    }
}

