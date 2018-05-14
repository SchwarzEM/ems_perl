#!/usr/bin/env perl

use strict;
use warnings;

my %term2number = ( 'ACEY.L3i'               =>  1,
                    'ACEY.24HCM'             =>  2,
                    'ACEY.24.PI'             =>  3,
                    'ACEY.5.D'               =>  4,
                    'ACEY.12.D'              =>  5,
                    'ACEY.17.D'              =>  6,
                    'ACEY.18D.Alb.4hr'       =>  7,
                    'ACEY.18D.noAlb.4hr'     =>  8,
                    'ACEY.18D.Cry5B.4hr'     =>  9,
                    'ACEY.18D.HEPES.4hr'     => 10,
                    'ACEY.18D.SB.plusCry5B'  => 11,
                    'ACEY.18D.SB.plusHEPES'  => 12,
                    'ACEY.18D.Cry5B.24hr'    => 13,
                    'ACEY.18D.HEPES.24hr'    => 14,
                    'ACEY.19.D'              => 15, );

print 'library(edgeR)', "\n";
print 'setwd("/sternlab/redivivus/data02/schwarz/Acey_genomics/post_meltdown/edgeR_2013.03.09.01")', "\n";
print 'Acey_hkeep_counts         <- read.delim("Acey_pme_expected_count.rounded.5plus_reads.01mar2013.406genes.txt",row.names="Gene")', "\n";
print 'Acey_hkeep_counts.DGEList <- DGEList(counts=Acey_hkeep_counts)', "\n";
print 'Acey_hkeep_counts.DGEList <- calcNormFactors(Acey_hkeep_counts.DGEList)', "\n";
print 'Acey_hkeep_counts.DGEList <- estimateCommonDisp(Acey_hkeep_counts.DGEList,verbose=TRUE)', "\n";

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

        my @initial_column_group_vals = (1..15);
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
        my $file_prefix = 'Acey_counts.' .  $control_term_text . q{.vs.} . $exptal_term;
        print "$file_prefix.counts",   ' <- read.delim("Acey_pme_expected_count.rounded.5plus_reads.01mar2013.txt",row.names="Gene")', "\n";
        print "$file_prefix.grouping", ' <- factor(c(', $final_col_grp_val_txt, '))', "\n";
        print "$file_prefix.DGEList",  ' <- DGEList(counts=', "$file_prefix.counts", ',group=', "$file_prefix.grouping", ") \n";
        print "$file_prefix.DGEList",  ' <- calcNormFactors(', "$file_prefix.DGEList", ") \n";
        print "$file_prefix.DGEList",  '$common.dispersion <- Acey_hkeep_counts.DGEList$common.dispersion', "\n";
        print "$file_prefix.DGEList.exactTest", ' <- exactTest(', "$file_prefix.DGEList, pair=c(1,2))\n";
        print "$file_prefix.DGEList.exactTest.topTags", ' <- topTags(', "$file_prefix.DGEList.exactTest, n=Inf)\n";
        print 'write.table(', "$file_prefix", '.DGEList.exactTest.topTags$table, file="', "$file_prefix", '.DGEList.exactTest.topTags.txt")', "\n";
    }
}

print "q()\n";

