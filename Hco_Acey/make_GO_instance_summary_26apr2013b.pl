#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my @ontologies = qw( molecular_function cellular_component biological_process );

my @data_sets = qw( acey_dirs cel_b2go_dirs hco_Pasi_dirs );

my %data_synonym = ( 'acey_dirs'     => "A. ceylanicum", 
                     'cel_b2go_dirs' => "C. elegans, blast2go",
                     'cel_wbGO_dirs' => "C. elegans, WS230 GO terms",
                     'hco_Pasi_dirs' => "H. contortus, gene predictions",
                     'hco_aug_dirs'  => "H. contortus, AUGUSTUS gene predictions", );

my %data_prefix = ( 'acey_dirs'     => "Acey_",
                    'cel_b2go_dirs' => "Cel_",
                    'cel_wbGO_dirs' => "Cel_GO.WS230_",
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
);

@{ $data_ref->{'cel_b2go_dirs'} } = qw ( 
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17feb2013.b2go.log10.L2.vs.L1.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.L3.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.herm_L4.vs.L3.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.YA.vs.herm_L4.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_17mar2013.b2go.log10.dauer_entry.vs.L2.func_wilcoxon_outdir
    /sternlab/redivivus/data02/schwarz/Acey_genomics/rsem_31oct2012/data_subsets/Cel_RSEM_11dec2012.b2go.log10.dauer_exit.vs.dauer_entry.func_wilcoxon_outdir
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
        my $prefix = $data_prefix{$data_set};
        # Sample $dir_basename:
        # Acey_RSEM_03dec2012.log10.24.PI.vs.L3i.func_wilcoxon_outdir
        if ( $dir_basename =~ /\A \S+ _RSEM \S+ \.log10\. (\S+) \.vs\. (\S+) \.func_wilcoxon_outdir \z/xms ) { 
            my $numerator   = $1;
            my $denominator = $2;
            $numerator   = $prefix . $numerator;
            $denominator = $prefix . $denominator;

            my $upreg_gene_command = 'filter_func_refinement.pl'
                                     . " -i $data_dir"
                                     . '/refinement-*-0.001_0.001.txt -u | cut -f 1-3,8'
                                     ;
            my $downreg_gene_command = 'filter_func_refinement.pl'    
                                     . " -i $data_dir"
                                     . '/refinement-*-0.001_0.001.txt -d | cut -f 1-3,7'
                                     ;
            my @upreg_term_text    = grep { $_ !~ /p_high_ranks_after_refinement/ } `$upreg_gene_command`;
            my @downreg_term_text  = grep { $_ !~ /p_low_ranks_after_refinement/ } `$downreg_gene_command`;
            chomp @upreg_term_text;
            chomp @downreg_term_text;

            foreach my $upreg_term (@upreg_term_text) { 
                if ( $upreg_term =~ /\A ([^\t]+) \t ([^\t]+ \t [^\t]+) \t (\S+) \z/xms ) {
                    my $ontology = $1;
                    my $go_term  = $2;
                    my $p_value  = $3;
                    my $instance = "$denominator to $numerator [+:";
                    $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'instance'}->{$instance}->{'p-value'} = $p_value;
                    if ( exists $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} ) { 
                        if ( $p_value < $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} ) {
                            $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} = $p_value;
                        }
                    }
                    if (! exists $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} ) { 
                        $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} = $p_value;
                    }
                }
                else {
                    die "Can't parse upregulated term $upreg_term\n";
                }
            }

            foreach my $downreg_term (@downreg_term_text) {
                if ( $downreg_term =~ /\A ([^\t]+) \t ([^\t]+ \t [^\t]+) \t (\S+) \z/xms ) {
                    my $ontology = $1;
                    my $go_term  = $2;
                    my $p_value  = $3;
                    my $instance = "$denominator to $numerator [-:";
                    $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'instance'}->{$instance}->{'p-value'} = $p_value;
                    if ( exists $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} ) {
                        if ( $p_value < $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} ) {
                            $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} = $p_value;
                        }
                    }
                    if (! exists $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} ) {
                        $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'min_p-value'} = $p_value;
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

foreach my $ontology (@ontologies) { 
    my @go_terms = sort {     $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$a}->{'min_p-value'} 
                          <=> $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$b}->{'min_p-value'} }
                   keys %{ $data_ref->{'ontology'}->{$ontology}->{'go_term'} };
    foreach my $go_term (@go_terms) { 
        my $markup = q{};
        my @instances = sort {     $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'instance'}->{$a}->{'p-value'}
                               <=> $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'instance'}->{$b}->{'p-value'} } 
                        keys %{ $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'instance'} };
        @instances = map { "$_ $data_ref->{'ontology'}->{$ontology}->{'go_term'}->{$go_term}->{'instance'}->{$_}->{'p-value'}]" } @instances;
        my $instance_text = join '; ', @instances;
        if (     ( $instance_text =~ / Acey_L3i \s+ to \s+ Acey_24HCM \s+ \[ [-] /xms ) 
             and ( $instance_text =~ / Cel_L3 \s+ to \s+ Cel_herm_L4 \s+ \[ [-] /xms )
             and ( $instance_text =~ / Hco_L3 \s+ to \s+ Hco_L4\.F \s+ \[ [-] /xms )
            ) {
            $markup = '---';
        }
        print "$ontology\t$markup\t$go_term\t$instance_text\n";
    }
}

