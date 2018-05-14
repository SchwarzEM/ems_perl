#!/usr/bin/env perl

# wilcoxon_gene_annots.pl -- Erich Schwarz <ems394@cornell.edu>, 4/17/2014.
# Purpose: given 1+ files with gene rankings and 1+ files with gene features, compute rank-sum (Wilcoxon) statistics for each feature.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(uniq);
use List::Util qw(max);
use Statistics::Descriptive;
use Statistics::Test::WilcoxonRankSum;

my @rankings  = ();
my @annots    = ();
my $qval_prog = q{};
my $data_ref;
my $help;

GetOptions ( 'rankings=s{,}' => \@rankings,
             'annots=s{,}'   => \@annots,
             'qval_prog=s'   => \$qval_prog,
             'help'          => \$help, );

if ($help or (! @rankings) or (! @annots)) { 
    die "Format: wilcoxon_gene_annots.pl\n",
        "            --rankings|-r    [1+ files or streams with gene rankings; e.g., log10(expr1/expr2)]\n",
        "            --annots|-a      [1+ files or streams with gene annots; 2+ annots/gene*column are split by default '; ']\n",
        '            --qval_prog|-q   [tell program where to find MEME\'s qvalue utility; if in $PATH, \'qvalue\' will do; else, give full path]', "\n",
        "            --help|-h        [print this message]\n",
        ;
}

# Standard format assumed for rankings files is like this:
#
# Gene	log10(24.PI/L3i)	log10(24.PI/24HCM)	log10(L3i/24.PI)
# Acey_2012.08.05_0001.g1	0.6112	0.4879	-0.6112
# Acey_2012.08.05_0001.g10	-0.2218	-0.4667	0.2218
# [etc.]

# Standard format assumed for annots files is like this:
#
# Gene	Pfam-A	Pfam-B	InterPro
# Acey_2012.08.05_0001.g131	3Beta_HSD [PF01073.14]; Epimerase [PF01370.16]; NAD_binding_10 [PF13460.1]; NmrA [PF05368.8]; RmlD_sub_bind [PF04321.12]		NAD(P)-binding domain [IPR016040]
# Acey_2012.08.05_0001.g122               Pfam-B_606 [PB000606]; Pfam-B_7725 [PB007725]
# Acey_2012.08.05_0002.g1087      SnoaL_3 [PF13474.1]     Pfam-B_15675 [PB015675]
# Acey_2012.08.05_0001.g100
# [etc.]

foreach my $ranking_file (@rankings) {
    # Accept either a stream from '-' or a standard file.

    my $INPUT_FILE;
    my $header = q{};
    my @rank_types = ();

    if ($ranking_file eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $ranking_file or die "Can't open input ranking file $ranking_file. $!\n";
    }

    # Record the incoming FASTA data.
    while (my $input_line = <$INPUT_FILE>) {
        chomp $input_line;

        # For each file, the very first input line is crucial because it tells us what to call the different rankings.
        if (! $header) { 
            $header = $input_line;
            if ( $header =~ /\A Gene ( (?: \t [^\t]+)+ ) \z/xms ) { 
                my $rank_text = $1;
                $rank_text =~ s/\A\t//;
                @rank_types = split /\t/, $rank_text;
            }
            else { 
                die "From ranking file $ranking_file, can't parse header line: $header\n";
            }
        }

        # Once we've read the first line that way, all subsequent lines should be giving us numbers that we can use to rank genes.
        else { 
            if ( $input_line =~ /\A ([^\t]+) ( (?: \t [^\t]+)+ ) \z/xms ) {
                my $gene = $1;
                my $rank_text = $2;
                $rank_text =~ s/\A\t//;
                my @rank_values = split /\t/, $rank_text;
                my $rank_index = @rank_values;
                $rank_index = $rank_index - 1;
                foreach my $rank_i (0..$rank_index) {
                    my $rank_type  = $rank_types[$rank_i];
                    my $rank_value = $rank_values[$rank_i];
                    if (! looks_like_number($rank_value) ) { 
                        die "Non-numeric rank \"$rank_value\" in data line of ranking file $ranking_file: $input_line\n";
                    }
                    $data_ref->{'rank_type'}->{$rank_type}->{'gene'}->{$gene}->{'rank_value'} = $rank_value;
                }
            }
            else {
                die "From ranking file $ranking_file, can't parse data line: $input_line\n";
            }
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to ranking file $ranking_file. $!\n";
}

foreach my $annot_file (@annots) {
    # Accept either a stream from '-' or a standard file.
    my $INPUT_FILE;
    my $header = q{};
    if ($annot_file eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }    
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $annot_file or die "Can't open input annotation file $annot_file. $!\n";
    }         

    # Record the incoming FASTA data.
    while (my $input_line = <$INPUT_FILE>) {
        chomp $input_line;

        # For each input file/etc., *ignore* first line, since we are not actually trying to remember what annnots came out of what column!
        if (! $header) { 
            $header = $input_line;
        }

        # After that very first line, we will assign individual annotations to genes willy-nilly.
        else { 
            # Crucial difference: it is OK, and in fact quite likely, that a given field will be devoid of annotations (of a given type) for any given gene.
            # So, "(\t [^\t]*)+", not "(\t [^\t]+)+".
            if ( $input_line =~ /\A ([^\t]+) ( (?: \t [^\t]*)+ ) \z/xms ) {
                my $gene = $1;
                my $initial_annot_text = $2;
                $initial_annot_text =~ s/\A\t//;

                my @initial_annots = split /\t/, $initial_annot_text;
                my @final_annots   = ();
                foreach my $initial_annot (@initial_annots) {
                    my @revised_annots = split /; /, $initial_annot;
                    push @final_annots, @revised_annots;
                }

                # If we somehow got 2+ copies of an annot, the next step will make them non-redundant again.
                # 
                # At the same time, for each annot, we also want to list the genes that have it,
                # so that we can eventually return that number to the statistics report.
                foreach my $final_annot (@final_annots) {
                    $data_ref->{'gene'}->{$gene}->{'annot'}->{$final_annot} = 1;
                    $data_ref->{'annot'}->{$final_annot}->{'gene'}->{$gene} = 1;
                }
            }
            else {
                die "From annotation file $annot_file, can't parse data line: $input_line\n";
            }
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to annotation file $annot_file. $!\n";
}

my @rank_types      = sort keys %{ $data_ref->{'rank_type'} };
my $rank_type_count = @rank_types;

foreach my $rank_type (@rank_types) { 
    # Sort genes by increasing rank values:
    my @genes = sort {     $data_ref->{'rank_type'}->{$rank_type}->{'gene'}->{$a}->{'rank_value'}
                       <=> $data_ref->{'rank_type'}->{$rank_type}->{'gene'}->{$b}->{'rank_value'} } 
                keys %{ $data_ref->{'rank_type'}->{$rank_type}->{'gene'} };

    my $gene_count = @genes;

    my @relevant_annots = ();
    foreach my $gene (@genes) {
        my @gene_annots = keys %{ $data_ref->{'gene'}->{$gene}->{'annot'} };
        push @relevant_annots, @gene_annots;
    }
    @relevant_annots = sort @relevant_annots;
    @relevant_annots = uniq @relevant_annots;

    my $annot_count = @relevant_annots;

    my @results   = ();
    foreach my $relevant_annot (@relevant_annots) {
        # Get the number of genes that had that annotation; we will print this later, as part of the overall statistics 
        #    (so that it is easy to tell if a given p-/q-value, for a given annotation, arose from 100 genes or 10).
        my $annot_gene_count = keys %{ $data_ref->{'annot'}->{$relevant_annot}->{'gene'} };

        # Set up to do the Wilcoxon rank-sum statistics.
        my @dataset_1 = ();
        my @dataset_2 = ();

        foreach my $gene (@genes) { 
            my $rank_value = $data_ref->{'rank_type'}->{$rank_type}->{'gene'}->{$gene}->{'rank_value'};
            if ( exists $data_ref->{'gene'}->{$gene}->{'annot'}->{$relevant_annot} ) {
                push @dataset_2, $rank_value;
            }
            else { 
                push @dataset_1, $rank_value;
            }
        }

        # The rank-sum test, by itself, does not tell us whether a strong p-value comes from *up* or *down*regulation; important to know this!
        # So, get an efficient average sign of change by taking the mean of log10(ratio), then reporting it for each group being assessed statistically.
        my $stat_dataset_2 = Statistics::Descriptive::Full->new();
        $stat_dataset_2->add_data(@dataset_2);
        my $mean_dataset_2 = $stat_dataset_2->mean();

        # Interpose this filter to avoid having the program choke on all-zero data sets.
        my $max_dataset_1 = max @dataset_1;
        my $max_dataset_2 = max @dataset_2;
        if ( ( $max_dataset_1 > 0 ) and ( $max_dataset_2 > 0 ) ) { 
            # Cribbed straight from perldoc Statistics::Test::WilcoxonRankSum:
            my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
            $wilcox_test->load_data(\@dataset_1, \@dataset_2);
            my $prob = $wilcox_test->probability();
            my @output = ("$rank_type", "$relevant_annot", "$annot_gene_count", "$mean_dataset_2", "$prob");
            push @results, \@output;
        }
    }

    # Generate a p-value to q-value mapping for this particular data set.
    # Doing so requires the qvalue utility program from MEME (http://meme.nbcr.net/meme/doc/qvalue.html).
    if ($qval_prog) {
        my @p_vals = map { $_->[4] } @results;

        # Purge '(', ')', and '/' [!] from a label that will go into a filename.
        my $rank_type_label = $rank_type;
        $rank_type_label    =~ s/[(]/./g;
        $rank_type_label    =~ s/[)]/./g;
        $rank_type_label    =~ s/\//.vs./g;

        my $tmp_pval_file = $rank_type_label . '.p_val_list.tmp.txt';
        my $tmp_qval_file = $rank_type_label . '.q_val_list.tmp.txt';

        $tmp_pval_file    =~ s/[.]+/./g;
        $tmp_qval_file    =~ s/[.]+/./g;

        $tmp_pval_file    = safename($tmp_pval_file);
        $tmp_qval_file    = safename($tmp_qval_file);

        open my $TMP_PVALS, '>', $tmp_pval_file or die "Can't open temp. p-value file \"$tmp_pval_file\"\n";
        foreach my $p_val (@p_vals) {
            print $TMP_PVALS "$p_val\n";
        }
        close $TMP_PVALS or die "Can't close filehandle to temp. p-value file $tmp_pval_file\n";

        system "$qval_prog --verbosity 1 $tmp_pval_file > $tmp_qval_file";

        open my $TMP_QVALS, '<', $tmp_qval_file or die "Can't open temp. q-value file \"$tmp_qval_file\"\n";
        while (my $input = <$TMP_QVALS>) { 
            chomp $input;
            if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
                my $p_val_x = $1;
                my $q_val_x = $2;
                $data_ref->{'rank_type'}->{$rank_type}->{'p_value'}->{$p_val_x}->{'q_value'} = $q_val_x;
            }
            else {
                die "Can't parse p- to q-value table $tmp_qval_file, at line: $input\n";
            }
        }
        close $TMP_QVALS or die "Can't close filehandle to temp. q-value file \"$tmp_qval_file\"\n";
    }

    @results = sort { $a->[4] <=> $b->[4] } @results; 
    foreach my $result (@results) {
        my $p_val_result = $result->[4];
        my $q_val_result = q{};
        if ( exists $data_ref->{'rank_type'}->{$rank_type}->{'p_value'}->{$p_val_result}->{'q_value'} ) { 
            $q_val_result = $data_ref->{'rank_type'}->{$rank_type}->{'p_value'}->{$p_val_result}->{'q_value'};
        }
        my $output = join "\t", @{ $result };
        $output = $output . "\t$q_val_result";
        print "$output\n";
    }
    $rank_type_count--;
    if ( $rank_type_count > 0 ) {
        print "\n";
    }
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

