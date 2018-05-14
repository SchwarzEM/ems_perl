#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;
use List::MoreUtils qw(uniq);

my $data_ref;

my @input_files = @ARGV;

foreach my $full_infile (@input_files) { 
    my $infile = basename $full_infile;
    if ( $infile =~ /\A AceyNoiseq\.(ACEY\S+)\.to\.(ACEY\S+)\.deg\.(?:up|down)\.raw\.txt \z/xms ) { 
        my $cond_a = $1;
        my $cond_b = $2;
        my $transition = "$cond_a to $cond_b";
        open my $IN, '<', $full_infile;
        while (my $input = <$IN>) {
            chomp $input;
            # Sample inputs:
            #   "ACEY.24.PI_mean" "ACEY.L3i_mean" "M" "D" "prob" "ranking"
            # "Acey_s0001.g130" 6252136.12558703 666.097749490496 13.1963276619642 6251470.02783754 1 6251470.02785147
            if ( ( $input !~ /ranking/ ) and ( $input =~ /\A \" (\S+) \" (?: \s+ \S+){2} \s+ (\S+) (?: \s+ \S+) \s+ (\S+) \s+ \S+ \s* \z/xms ) ) { 
                my $gene  = $1;
                my $fold  = $2;
                my $p_val = $3;

                my $sign  = q{};

                if ($fold > 0) {
                    $sign = '+';
                }
                if ($fold < 0) {
                    $sign = '-';
                }
                if ($fold == 0) {
                    die "Can't parse zero-fold log change of activity in input file $full_infile: $input\n";
                }

                # I want p-value sorted outputs like this, per gene: ACEY.5.D to ACEY.12.D [+: 0.000104908789422638]
                my $transition = "$cond_a to $cond_b" . q{ [} . $sign . ": $p_val" . q{]};

                # For each gene, we store all our annotations in processed form, accessible by their p-values.
                $data_ref->{'gene'}->{$gene}->{'p_value_data'}->{$p_val}->{'transition'}->{$transition} = 1;

                # But we want their p-values to be easily sortable numerically, which requires a *separate* reference.
                push @{ $data_ref->{'gene'}->{$gene}->{'p_values'} }, $p_val;
            }
        }
    }
    else {
        die "Can't parse input file: $full_infile\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };

my $header = "Gene\tNOISeq";

foreach my $gene (@genes) { 
    my @annots = ();
    print "$header\n" if $header;
    $header = q{};

    # For NOISeq, the p-value is *better* if it is *higher* (closer to 1.0 than 0.0), 
    #     so we sort by $b <=> $a rather than the reverse (that we would use for DESeq q-values).
    my @p_values = sort { $b <=> $a } @{ $data_ref->{'gene'}->{$gene}->{'p_values'} };
    @p_values = uniq(@p_values);
    foreach my $p_value (@p_values) {
        my @new_annots = sort keys %{ $data_ref->{'gene'}->{$gene}->{'p_value_data'}->{$p_value}->{'transition'} };
        push @annots, @new_annots;
    }

    my $annot_text = join '; ', @annots;
    print "$gene\t$annot_text\n";
}

