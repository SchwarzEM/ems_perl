#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my @data_sets = qw( acey_dirs cel_b2go_dirs cel_wbGO_dirs hco_Pasi_dirs hco_aug_dirs );

my %data_synonym = ( 'acey_dirs'     => "A. ceylanicum", 
                     'cel_b2go_dirs' => "C. elegans, blast2go GO terms",
                     'cel_wbGO_dirs' => "C. elegans, WS230 GO terms",
                     'hco_Pasi_dirs' => "H. contortus, gene predictions",
                     'hco_aug_dirs'  => "H. contortus, AUGUSTUS gene predictions", );

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

            my $upreg_gene_command = 'filter_func_refinement.pl'
                                     . " -i $data_dir"
                                     . '/refinement-*-0.001_0.001.txt -u | cut -f 1-3,8'
                                     ;
            my $downreg_gene_command = 'filter_func_refinement.pl'    
                                     . " -i $data_dir"
                                     . '/refinement-*-0.001_0.001.txt -d | cut -f 1-3,7'
                                     ;
            my @upreg_term_text    = ();

            my @downreg_term_text    = ();

            my @go_ontologies     = qw( molecular_function cellular_component biological_process );
            foreach my $go_ont (@go_ontologies) { 
                my @upreg_scratch_text = sort { get_p_value($a) <=> get_p_value($b) }
                                         grep { $_ =~ /\A$go_ont/ } 
                                         grep { $_ !~ /p_high_ranks_after_refinement/ } 
                                         `$upreg_gene_command`;
                push @upreg_term_text, @upreg_scratch_text;
                push @upreg_term_text, "\n";
            }
            foreach my $go_ont (@go_ontologies) {
                my @downreg_scratch_text = sort { get_p_value($a) <=> get_p_value($b) }
                                           grep { $_ =~ /\A$go_ont/ } 
                                           grep { $_ !~ /p_low_ranks_after_refinement/ } 
                                           `$downreg_gene_command`;
                push @downreg_term_text, @downreg_scratch_text;
                push @downreg_term_text, "\n";
           }

            print "GO terms associated with upregulation from $denominator to $numerator in $data_synonym{$data_set}:\tXXXXX\tXXXXX\tXXXXX\nXXXXX\tXXXXX\tXXXXX\tXXXXX\n";
            print "root_node_name\tnode_name\tnode_id\tp_high_ranks_after_refinement\tXXXXX\tXXXXX\tXXXXX\n";
            print @upreg_term_text;
            print "XXXXX\tXXXXX\tXXXXX\tXXXXX\n";
            print "GO terms associated with downregulation from $denominator to $numerator in $data_synonym{$data_set}:\tXXXXX\tXXXXX\tXXXXX\nXXXXX\tXXXXX\tXXXXX\tXXXXX\n";
            print "root_node_name\tnode_name\tnode_id\tp_low_ranks_after_refinement\nXXXXX\tXXXXX\tXXXXX\tXXXXX\n";
            print @downreg_term_text;
            print "XXXXX\tXXXXX\tXXXXX\tXXXXX\n";
        }
        else {
            die "Can't parse specific label in $dir_basename\n";
        }
    }
}

sub get_p_value {
    my $_input_line = $_[0];
    chomp $_input_line;
    if ( $_input_line =~ /\A [^\t]+ \t [^\t]+ \t [^\t]+ \t (\S+) \z/xms ) { 
        my $_p_value = $1;
        return $_p_value;
    }
    else { 
        die "Failed to parse p-value of $_input_line\n";
    }
}

