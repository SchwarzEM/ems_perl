#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my $positive_raw_data = $ARGV[0];
my $negative_raw_data = $ARGV[1];

my $data_ref;

process_data($positive_raw_data, 'positive');
process_data($negative_raw_data, 'negative');

my @all_positive_genes = sort keys %{ $data_ref->{'class'}->{'positive'}->{'gene'} };
my @chosen_pos_genes   = ();

foreach my $all_pos_gene (@all_positive_genes) {
    if (! $data_ref->{'class'}->{'positive'}->{'gene'}->{$all_pos_gene}->{'expr_val'} ) {
        die "Cannot identify positive expression value for $all_pos_gene\n";
    }
    my $pos_gene_expr_val = $data_ref->{'class'}->{'positive'}->{'gene'}->{$all_pos_gene}->{'expr_val'};

    if (! exists $data_ref->{'class'}->{'negative'}->{'gene'}->{$all_pos_gene}->{'expr_val'} ) {
        push @chosen_pos_genes, $all_pos_gene;
    }
    else { 
        my $neg_gene_expr_val = $data_ref->{'class'}->{'negative'}->{'gene'}->{$all_pos_gene}->{'expr_val'};
        if ( $neg_gene_expr_val == 0 ) {
            push @chosen_pos_genes, $all_pos_gene;
        }
        my $gene_exp_ratio = ( $pos_gene_expr_val / $neg_gene_expr_val );
        if ( $gene_exp_ratio >= 3 ) {
            push @chosen_pos_genes, $all_pos_gene;
        }
    }
}

@chosen_pos_genes = sort @chosen_pos_genes;
@chosen_pos_genes = uniq @chosen_pos_genes;

foreach my $chosen_pos_gene (@chosen_pos_genes) {
    my $pos_gene_expr_val = $data_ref->{'class'}->{'positive'}->{'gene'}->{$chosen_pos_gene}->{'expr_val'};
    print "$chosen_pos_gene\t$pos_gene_expr_val\n";
}


sub process_data {
    my $_raw_data_file  = $_[0];
    my $_raw_data_class = $_[1];
    open my $_RAW, '<', $_raw_data_file;
    while (my $_input = <$_RAW>) {
        chomp $_input;
        if ( $_input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+\S*) \b .* \t (\S+) \z/xms ) { 
            my $_gene_text = $1;
            my $_expr_val  = $2;

            if (! looks_like_number($_expr_val) ) {
                die "Non-numerical putative expression level value in $_raw_data_file for: $_input\n";
            }

            # because these dork files sometimes have 'AT1G07590;AT1G07600' names
            my @_genes = split /;/, $_gene_text;

            foreach my $_gene (@_genes) {
                # These goofy files actually *do* have redundant data!
                # Pragmatic solution: pick the largest expression value for any given gene that the data have.
                if ( exists $data_ref->{'class'}->{$_raw_data_class}->{'gene'}->{$_gene}->{'expr_val'} ) {
                    my $_prev_expr_val = $data_ref->{'class'}->{$_raw_data_class}->{'gene'}->{$_gene}->{'expr_val'} ;
                    if ( $_prev_expr_val > $_expr_val ) {
                        $_expr_val = $_prev_expr_val;
                    }
                }
                $data_ref->{'class'}->{$_raw_data_class}->{'gene'}->{$_gene}->{'expr_val'} = $_expr_val;
            }
        }
    }
    close $_RAW;
    return;
}
