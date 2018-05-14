#!/usr/bin/env perl

use strict;
use warnings;

my %term2number = ( 'Cel_early_embryo'        => 1,
                    'Cel_embryo_0min'         => 2,
                    'Cel_embryo_30min'        => 3,
                    'Cel_embryo_60min'        => 4,
                    'Cel_embryo_90min'        => 5,
                    'Cel_embryo_120min'       => 6,
                    'Cel_late_embryo'         => 7,
                    'Cel_L1'                  => 8,
                    'Cel_L2'                  => 9, 
                    'Cel_L3'                  => 10,
                    'Cel_herm_L4'             => 11,
                    'Cel_Alb.4hr.L4'          => 12,
                    'Cel_noAlb.cont.L4'       => 13,
                    'Cel_male_L4'             => 14,
                    'Cel_dauer_entry'         => 15,
                    'Cel_dauer_exit'          => 16,
                    'Cel_YA'                  => 17, );

print 'library(edgeR)', "\n";
print 'setwd("/sternlab/redivivus/data02/schwarz/Acey_genomics/post_meltdown/edgeR_2013.03.09.01")', "\n";
print 'Cel_hkeep_counts         <- read.delim("Cel_modencode_pme_exp_count.5plus_reads.rounded.01mar2013.406genes.txt",row.names="Gene")', "\n";
print 'Cel_hkeep_counts.DGEList <- DGEList(counts=Cel_hkeep_counts)', "\n";
print 'Cel_hkeep_counts.DGEList <- calcNormFactors(Cel_hkeep_counts.DGEList)', "\n";
print 'Cel_hkeep_counts.DGEList <- estimateCommonDisp(Cel_hkeep_counts.DGEList,verbose=TRUE)', "\n";

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    if ( $input =~ /\A (\S.+\S) \s+ vs\. \s+ (\S.+\S) \z/xms ) { 
        my $control_text   = $1;
        my $exptal_term    = $2;
        my %is_control_val = ();

        my @control_terms = split /\s+[+]\s+/, $control_text;
        my $control_term_text = join '_and_', @control_terms;
        my @control_vals  = map { $term2number{$_} } @control_terms;
        foreach my $control_val (@control_vals) { 
            $is_control_val{$control_val} = 1;
        }
        my $exptal_val    = $term2number{$exptal_term};

        my @initial_column_group_vals = (1..17);
        my @final_column_group_vals    = ();
        foreach my $init_col_grp_val (@initial_column_group_vals) { 
            if ( $init_col_grp_val == $exptal_val ) { 
                push @final_column_group_vals, '2';
            }
            elsif ( exists $is_control_val{$init_col_grp_val} ) { 
                push @final_column_group_vals, '1';
            }
            else { 
                $init_col_grp_val += 20;
                push @final_column_group_vals, $init_col_grp_val;
            }
        }       
        my $final_col_grp_val_txt = join q{,}, @final_column_group_vals;
        my $file_prefix = 'Cel_counts.' .  $control_term_text . q{.vs.} . $exptal_term;
        print "$file_prefix.counts",   ' <- read.delim("Cel_modencode_pme_exp_count.5plus_reads.rounded.01mar2013.txt",row.names="Gene")', "\n";
        print "$file_prefix.grouping", ' <- factor(c(', $final_col_grp_val_txt, '))', "\n";
        print "$file_prefix.DGEList",  ' <- DGEList(counts=', "$file_prefix.counts", ',group=', "$file_prefix.grouping", ") \n";
        print "$file_prefix.DGEList",  ' <- calcNormFactors(', "$file_prefix.DGEList", ") \n";
        print "$file_prefix.DGEList",  '$common.dispersion <- Cel_hkeep_counts.DGEList$common.dispersion', "\n";
        print "$file_prefix.DGEList.exactTest", ' <- exactTest(', "$file_prefix.DGEList, pair=c(1,2))\n";
        print "$file_prefix.DGEList.exactTest.topTags", ' <- topTags(', "$file_prefix.DGEList.exactTest, n=Inf)\n";
        print 'write.table(', "$file_prefix", '.DGEList.exactTest.topTags$table, file="', "$file_prefix", '.DGEList.exactTest.topTags.txt")', "\n";
    }
}

print "q()\n";

