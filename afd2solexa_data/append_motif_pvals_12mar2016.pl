#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $p_val_data = q{};
my $motif_data = q{};

my $data_set   = q{};
my $motif      = q{};
my $p_val_type = q{};
my $p_val      = q{};

my $added_head = "Query_p-value\tSepal_p-value";

if ($ARGV[0] and $ARGV[1]) {
    $p_val_data = $ARGV[0];
    $motif_data = $ARGV[1];
}
else {
    die "Format: append_motif_pvals_12mar2016.pl [p-value data] [motif data]\n";
}

open my $PVAL_DAT, '<', $p_val_data;
while (my $input = <$PVAL_DAT>) {
    chomp $input;
    # > # Data set: "all_sepal_geno_diff_expr_gene_list"; motif "1"; p-value for freq. in all_sepal_genes
    # > # Data set: "all_sepal_geno_diff_expr_gene_list"; motif "1"; p-value for freq. in all_query_genes
    if ( $input =~ /\A [>] [ ] [#] [ ] Data [ ] set: [ ] \" ([^\"]+) \" [;] 
                       [ ] motif [ ] \" ([^\"]+) \" [;] 
                       [ ] p-value [ ] for [ ] freq\. [ ] in [ ] (\S+) \s* \z/xms ) {
        $data_set   = $1;
        $motif      = $2;
        $p_val_type = $3;
    }
    elsif ( $data_set and ( $input =~ /\A number [ ] of [ ] successes [ ] = [ ]  \d+, 
                            [ ] number [ ] of [ ] trials [ ] = [ ]  \d+, 
                            [ ] p-value [ ] (?: <|=) [ ]  (\S+) \s* \z/xms ) ) {
        $p_val  = $1;
        $data_ref->{'data_set'}->{$data_set}->{'motif'}->{$motif}->{'p_val_type'}->{$p_val_type} = $p_val;
        $data_set   = q{};
        $motif      = q{};
        $p_val_type = q{};
        $p_val      = q{};
    }
}
close $PVAL_DAT;

open my $MOT_DAT, '<', $motif_data;
while (my $input = <$MOT_DAT>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) .* \z/xms ) {
        my $data_set1 = $1;
        my $motif1    = $2;

        if ( $data_set1 eq 'Data' ) {
            print "$input\t$added_head\n";
        }
        elsif (! exists $data_ref->{'data_set'}->{$data_set1}->{'motif'}->{$motif1}->{'p_val_type'} ) {
            die "Do not have p-value data for motif in: $input\n";
        }
        else {
            my $query_p_value = $data_ref->{'data_set'}->{$data_set1}->{'motif'}->{$motif1}->{'p_val_type'}->{'all_query_genes'};
            my $sepal_p_value = $data_ref->{'data_set'}->{$data_set1}->{'motif'}->{$motif1}->{'p_val_type'}->{'all_sepal_genes'};
            print "$input\t$query_p_value\t$sepal_p_value\n";
        }
    }
    else {
        die "From motif data table $motif_data, cannot parse: $input\n";
    }
}
close $MOT_DAT;

