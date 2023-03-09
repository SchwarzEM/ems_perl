#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);
use Text::NSP::Measures::2D::Fisher2::twotailed;

my $data_ref;

my $total_gene_count = 0;
my $class_desc       = q{};

my $header = "Motif\t"
             . "All_genes\t"
             . "Motif_genes\t"
             . "Class_genes\t"
             . "Motif.Class_overlap\t"
             . "Exp_rand_overlap\t"
             . "Enrichment\t"
             . "p-value\t"
             . "q-value"
             ;

my $data = q{};
$data    = $ARGV[0] if $ARGV[0];

my $qval_prog = q{};
$qval_prog    = $ARGV[1] if $ARGV[1];

my @results = ();

if ( (! -e $data ) or (! -e $qval_prog ) ) {
    die "Format: motif_group_fisher_07mar2023.pl",
        " [input TSV table: gene / motifs / Boolean class] [qvalue program]",
        " > [Fisher two-tailed p-values and q-values for motif/class overlaps]\n",
        ;
}

open my $DATA, '<', $data;
while (my $input = <$DATA>) {
    chomp $input;
    if ( ( $input =~ /\A (\S+) \t ([^\t]*) \t ([^\t]*) \z/xms ) and ( $input !~ /\AGene/xms ) ) {
        my $gene    = $1;
        my $mot_txt = $2;
        my $class   = $3;

        $total_gene_count++;

        if ( $mot_txt =~ /\S/xms ) {
            my @motifs    = split '; ', $mot_txt;
            @motifs    = sort @motifs;
            @motifs    = uniq @motifs;

            foreach my $motif (@motifs) {
                $data_ref->{'motif'}->{$motif}->{'gene'}->{$gene} = 1;
            }
        }

        if ( $class =~ /\S/xms ) {
            if ( $class_desc and ( $class_desc ne $class ) ) {
                die "Two different Boolean class descriptions: \"$class_desc\" versus \"$class\"\n";
            }
            $class_desc = $class;
            $data_ref->{'gene'}->{$gene}->{'class'}->{$class} = 1;
            $data_ref->{'class_gene'}->{$gene} = 1;
        }
    }
    elsif ( $input !~ /\AGene/xms ) {
        die "From data file $data, cannot parse: $input\n";
    }
}
close $DATA;

my @motifs = sort keys %{ $data_ref->{'motif'} };

foreach my $motif (@motifs) {
    my $total_class_genes = q{};
    $total_class_genes    = keys %{ $data_ref->{'class_gene'} };

    my $total_motif_genes   = q{};
    $total_motif_genes      = keys %{ $data_ref->{'motif'}->{$motif}->{'gene'} };

    my $motif_class_overlap = 0;
    my @genes = sort keys %{ $data_ref->{'motif'}->{$motif}->{'gene'} };
    foreach my $gene (@genes) {
        if ( exists $data_ref->{'gene'}->{$gene}->{'class'} ) {
            $motif_class_overlap++;
        }
    }
    my $n11 = $motif_class_overlap;
    my $n1p = ($total_motif_genes - $motif_class_overlap);
    my $np1 = ($total_class_genes - $motif_class_overlap);
    my $npp = $total_gene_count;

    my $fisher_pval = 'NA';
    my $f_score     = calculateStatistic( n11=>$n11,
                                          n1p=>$n1p,
                                          np1=>$np1,
                                          npp=>$npp);

    # Any motif that fails to get a Fisher p-value at all is discarded in this step:
    if ( looks_like_number($f_score) ) {
        $fisher_pval = $f_score;

        my $exp_rand_overlap = ($total_motif_genes * ($total_class_genes / $total_gene_count)); 
        my $enrichment       = ($motif_class_overlap / $exp_rand_overlap);

        my @output = ( "$motif", 
                       "$total_gene_count", 
                       "$total_motif_genes", 
                       "$total_class_genes", 
                       "$motif_class_overlap", 
                       "$exp_rand_overlap",
                       "$enrichment",
                       "$fisher_pval", 
                     );
        push @results, \@output;
    }
}

my @p_vals  = sort { $a <=> $b } map { $_->[7] } @results;
my @meme_ps = ();
my @q_vals  = ();

my $tmp_pval_file = $data . '.p_val_list.tmp.txt';
my $tmp_qval_file = $data . '.q_val_list.tmp.txt';
$tmp_pval_file    = safename($tmp_pval_file);
$tmp_qval_file    = safename($tmp_qval_file);

open my $TMP_PVALS, '>', $tmp_pval_file;
foreach my $p_val (@p_vals) {
    print $TMP_PVALS "$p_val\n";
}
close $TMP_PVALS;

# '--append' is crucial, because it passes on the original p-values unchanged to the p-/q-value table.
system "$qval_prog --append --verbosity 1 $tmp_pval_file > $tmp_qval_file";

open my $TMP_QVALS, '<', $tmp_qval_file;
while (my $input = <$TMP_QVALS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $p_val = $1;
        my $q_val = $2;
        $data_ref->{'p_value'}->{$p_val}->{'q_value'} = $q_val;
    }
    else {
        die "Can't parse p- to q-value table $tmp_qval_file, at line: $input\n";
    }
}
close $TMP_QVALS;

@results = sort { $a->[7] <=> $b->[7] } @results; 

foreach my $result (@results) {
        my $p_val  = $result->[7];
        my $q_val  = q{};

        if ( exists $data_ref->{'p_value'}->{$p_val}->{'q_value'} ) { 
            $q_val = $data_ref->{'p_value'}->{$p_val}->{'q_value'};
        }
        else {
            die "Failed to map p-value to q-value: $p_val\n";
        }

        my $output = join "\t", @{ $result };
        $output = $output . "\t$q_val";

        print "$header\n" if $header;
        $header = q{};

        print "$output\n";
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

