#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
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

my $data      = q{};
$data         = $ARGV[0] if $ARGV[0];

my $qval_prog = q{};
$qval_prog    = $ARGV[1] if $ARGV[1];

my $basename  = basename($data);

my @results   = ();

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

    my $fisher_pval = q{};
    $fisher_pval    = calculateStatistic( n11=>$motif_class_overlap,
                                          n1p=>$total_motif_genes,
                                          np1=>$total_class_genes,
                                          npp=>$total_gene_count, );

    my $errorCode    = q{};
    my $errorMessage = q{};
    if ( ($errorCode = getErrorCode()) ) {
         $errorMessage = getErrorMessage();
         warn "Motif $motif has error $errorCode: $errorMessage\n";
    }

    # Any motif that fails to get a Fisher p-value at all is discarded in this step:
    if ( looks_like_number($fisher_pval) ) {

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
    else {
        warn "Failed to assign Fisher p-value to motif: $motif\n";
        warn "Putative Fisher p-value of motif $motif was: $fisher_pval\n";
        warn "    Motif class overlap (n11) was: $motif_class_overlap\n";
        warn "    Total motif genes (n1p) was: $total_motif_genes\n";
        warn "    Total class genes (np1) was: $total_class_genes\n";
        warn "    Total gene count (npp) was: $total_gene_count\n";
    }
}

@results   = sort { $a->[7] <=> $b->[7] } @results;
my @p_vals = map { $_->[7] } @results;

my @meme_ps = ();
my @q_vals  = ();

my $tmp_pval_file = $basename . '.p_val_list.tmp.txt';
my $tmp_qval_file = $basename . '.q_val_list.tmp.txt';
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
    if ( $input =~ /\A \S+ \t (\S+) \z/xms ) {
        my $q_val = $1;
        push @q_vals, $q_val;
    }
    else {
        die "Can't parse p- to q-value table $tmp_qval_file, at line: $input\n";
    }
}
close $TMP_QVALS;

my $result_count = q{};
my $q_val_count  = q{};
$result_count    = @results;
$q_val_count     = @q_vals;

if ( $result_count != $q_val_count ) {
    die "Unequal numbers of results ($result_count) versus q-values ($q_val_count)\n";
}

$result_count--;

foreach my $i (0..$result_count) {
    my $result = $results[$i];
    my $q_val  = $q_vals[$i];

    my $output = join "\t", @{ $result };
    $output    = $output . "\t$q_val";

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

