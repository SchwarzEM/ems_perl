#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);
use Text::NSP::Measures::2D::Fisher2::twotailed;

my $data_ref;

my $data      = q{};
$data         = $ARGV[0] if $ARGV[0];

my $qval_prog = q{};
$qval_prog    = $ARGV[1] if $ARGV[1];

my $total_genes_A = q{};
$total_genes_A    = $ARGV[2] if $ARGV[2];

my $total_genes_B = q{};
$total_genes_B    = $ARGV[3] if $ARGV[3];

my $total_genes_both = q{};

if ( (! -e $data ) or (! -e $qval_prog ) or (! $total_genes_A ) or (! $total_genes_B ) ) {
    kill_program();
}

if ( (! looks_like_number($total_genes_A) ) or (! looks_like_number($total_genes_B) ) ) {
    warn "Either total gene count A \"$total_genes_A\" or total gene count B \"$total_genes_B\" was not a number.\n";
    kill_program();
}

$total_genes_both = $total_genes_A + $total_genes_B;

my $header = "Motif\t"
             . "Total_genes_A\t"
             . "Mot_genes_A\t"
             . "Total_genes_B\t"
             . "Mot_genes_B\t"
             . "Exp_rand_overlap\t"
             . "Enrichment\t"
             . "p-value\t"
             . "q-value"
             ;

my $basename  = basename($data);

my @results   = ();

open my $DATA, '<', $data;
while (my $input = <$DATA>) {
    chomp $input;
    if ( ( $input =~ /\A ([^\t]+) \t (\d+) \t (\d+) \z/xms ) and ( $input !~ /\AGene/xms ) ) {
        my $motif    = $1;
        my $mot_genes_A = $2;
        my $mot_genes_B = $3;
        my $mot_genes_both = $mot_genes_A + $mot_genes_B;

        if ( exists $data_ref->{'motif'}->{$motif}->{'mot_genes_A'} ) {
            die "Redundant counts of motif $motif: $data_ref->{'motif'}->{$motif}->{'mot_genes_A'} vs. $mot_genes_A\n";
        }
        if ( exists $data_ref->{'motif'}->{$motif}->{'mot_genes_B'} ) { 
            die "Redundant counts of motif $motif: $data_ref->{'motif'}->{$motif}->{'mot_genes_B'} vs. $mot_genes_B\n";
        }
        if ( exists $data_ref->{'motif'}->{$motif}->{'mot_genes_both'} ) {
            die "Redundant total count of motif $motif: $data_ref->{'motif'}->{$motif}->{'mot_genes_both'} vs. $mot_genes_both\n";
        }

        $data_ref->{'motif'}->{$motif}->{'mot_genes_A'}    = $mot_genes_A;
        $data_ref->{'motif'}->{$motif}->{'mot_genes_B'}    = $mot_genes_B;
        $data_ref->{'motif'}->{$motif}->{'mot_genes_both'} = $mot_genes_both;
    }
    elsif ( $input !~ /\AMotif/xms ) {
        die "From data file $data, cannot parse: $input\n";
    }
}
close $DATA;

my @motifs = sort keys %{ $data_ref->{'motif'} };

foreach my $motif (@motifs) {
    if (! exists $data_ref->{'motif'}->{$motif}->{'mot_genes_A'} ) {
        die "Cannot map motif $motif to mot_genes_A count\n";
    }
    if (! exists $data_ref->{'motif'}->{$motif}->{'mot_genes_B'} ) {
        die "Cannot map motif $motif to mot_genes_B\n";
    }
    if (! exists $data_ref->{'motif'}->{$motif}->{'mot_genes_both'} ) {
        die "Cannot map motif $motif to mot_genes_both\n";
    }

    my $mot_genes_A = $data_ref->{'motif'}->{$motif}->{'mot_genes_A'};
    my $mot_genes_B = $data_ref->{'motif'}->{$motif}->{'mot_genes_B'};
    my $mot_genes_both = $data_ref->{'motif'}->{$motif}->{'mot_genes_both'};

    my $fisher_pval = q{};
    $fisher_pval    = calculateStatistic( n11=>$mot_genes_A,
                                          n1p=>$mot_genes_both,
                                          np1=>$total_genes_A,
                                          npp=>$total_genes_both, );
    my $errorCode    = q{};
    my $errorMessage = q{};
    if ( ($errorCode = getErrorCode()) ) {
         $errorMessage = getErrorMessage();
         warn "Motif $motif has error $errorCode: $errorMessage\n";
    }


    # Any motif that fails to get a Fisher p-value at all is discarded in this step:
    if ( looks_like_number($fisher_pval) ) {

        my $exp_rand_overlap = ( $mot_genes_both * ($total_genes_A / $total_genes_both) ); 
        my $enrichment       = ($mot_genes_A / $exp_rand_overlap);

        my @output = ( "$motif", 
                       "$total_genes_A", 
                       "$mot_genes_A",
                       "$total_genes_B",
                       "$mot_genes_B",
                       "$exp_rand_overlap",
                       "$enrichment",
                       "$fisher_pval",
                     );
        push @results, \@output;
    }
    else {
        warn "Failed to assign Fisher p-value to motif: $motif\n";
        warn "Putative Fisher p-value of motif $motif was: $fisher_pval\n";
        warn "    Motif gene count in A (n11) was:   $mot_genes_A\n";
        warn "    Motif gene count in A+B (n1p) was: $mot_genes_both\n";
        warn "    Total gene count in A (np1) was:   $total_genes_A\n";
        warn "    Total gene count in A+B (npp) was: $total_genes_both\n";
    }
}

@results   = sort { $a->[7] <=> $b->[7] } @results;
my @p_vals = map { $_->[7] } @results;

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

sub kill_program {
    die "Format: motif_group_fisher_11sep2024.pl\n",
        "            [input TSV table: motif / motif-gene-count-A / motif-gene-count-B]\n",
        "            [qvalue program]\n",
        "            [total gene count, taxon A]\n",
        "            [total gene count, taxon B]\n",
        "        > [Fisher two-tailed p-values and q-values for motif/class overlaps]\n",
        ;
    return;
}

