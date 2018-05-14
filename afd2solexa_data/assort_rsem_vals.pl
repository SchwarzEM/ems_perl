#!/usr/bin/env perl

# assort_rsem_vals.pl -- Erich Schwarz <emsch@caltech.edu>, 11/26/2012.
# Purpose: given RSEM values in simple two-column files, stitch them into a single report with defined ratios and no irrelevant genes.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @reports_w_gene_data       = ();
my @labels_of_reports         = ();
my @reports_w_important_genes = ();
my $nozeros;
my @wanted_ratios             = ();
my $pseudozero                = 0;
my $help;

my $data_ref;

my $gene           = q{};
my $rsem           = q{};
my $i              = 0;

GetOptions ( 'data=s{,}'      => \@reports_w_gene_data,
             'labels=s{,}'    => \@labels_of_reports,
             'important=s{,}' => \@reports_w_important_genes,
             'nozeros'        => \$nozeros,
             'ratios=s{,}'    => \@wanted_ratios, 
             'pseudozero=f'   => \$pseudozero,
             'help'           => \$help,   );

my $data_count  = @reports_w_gene_data;

# Labels default to full filenames, if no label names are given:
if (! @labels_of_reports) {
    @labels_of_reports = @reports_w_gene_data;
}

# *After* we know we have something in @labels_of_reports, count it:
my $label_count = @labels_of_reports;

# Default is to consider genes for every data file 'important', and thus subject to reporting and non-zero filter:
if (! @reports_w_important_genes) {
    @reports_w_important_genes = @reports_w_gene_data;
}

if ( $help or (! @reports_w_gene_data) or ( $data_count != $label_count ) ) { 
    die "Format: assort_rsem_vals.pl\n",
        "    --data|-d         <N files, with RSEM data for genes>\n",
        "    --labels|-l       <N labels, to put atop columns of data from each listed file; if listed, must equal N; if not listed, defaults to data file names>\n",
        "    --important|-i    <labels of files whose genes we actually want to report; defaults to all N data files>\n",
        "    --nozeros|-n      <however the files listing reported genes were chosen, only list those with a nonzero RSEM value in that file set>\n",
        "    --ratios|-r       <ratios that we want (if any); list as [label1]/[label2]>\n",
        "    --pseudozero|-p   [set an arbitrary, universal, fixed pseudozero such as 0.01; the default is to the smallest observed nonzero for each data set]\n",
        "    --help|-h         [print this message]\n",
        ;
}

# Convert from 1-to-N to 0-to-(N-1) numbering.
$data_count--;

# This allows files to be identified by their labels (a lot less clunky than ID by full filenames!), or vice versa.
foreach my $i (0..$data_count) {
    my $label        = $labels_of_reports[$i];
    my $labeled_file = $reports_w_gene_data[$i];
    $data_ref->{'file'}->{$labeled_file}->{'file_abbrev'} = $label;
    $data_ref->{'file_abbrev'}->{$label}->{'file'} = $labeled_file;
}

foreach my $imp_gene_file (@reports_w_important_genes) {
    if (! exists $data_ref->{'file'}->{$imp_gene_file}->{'file_abbrev'} ) { 
        die "Can't identify which important label refers to important file $imp_gene_file\n";
    }
    open my $IMP_GENES, '<', $imp_gene_file or die "Can't open file $imp_gene_file with wanted gene names: $!";
    while (my $input = <$IMP_GENES>) { 
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
            $gene = $1;
            $rsem = $2;
            # Use looks_like_number to specifically avoid getting header lines or other irrelevancies entered as numerical data.
            if ( ( looks_like_number $rsem ) and ( (! $nozeros) or ( $rsem > 0 ) ) ) { 
                $data_ref->{'wanted_gene'}->{$gene} = 1;
            }
        }
    }
    close $IMP_GENES or die "Can't close filehandle to file $imp_gene_file with wanted gene names: $!";
}

# Having gotten all wanted gene names, get them into a list for later use:
my @genes_of_interest = sort keys %{ $data_ref->{'wanted_gene'} };

foreach my $i (0..$data_count) { 
    my $label = $labels_of_reports[$i];
    open my $DATA_FILE, '<', $reports_w_gene_data[$i] or die "Can't open data file $reports_w_gene_data[$i] (labeled $labels_of_reports[$i]): $!";
    while (my $input = <$DATA_FILE>) { 
        chomp $input;
        if ( $input =~ / \A (\S+) \t (\S+) \s* \z /xms) { 
            $gene = $1;
            $rsem = $2;

            # Sanity-check numbers, for genes which matter (i.e., we do not obsess over non-numerical lines or irrelevant genes).
            # Note that this error should not be possible, given previous use of looks_like_number as an input filer; so, argubly overdefensive coding here.
            if ( ( exists $data_ref->{'wanted_gene'}->{$gene} ) and (! ( looks_like_number $rsem ) ) ) {
                die "In $input: for gene $gene, RSEM value $rsem does not look numerical.\n";
            }

            # The business about overwriting lower values with larger ones is a relic from parsing funky WBGene retro-mapped data onto split gene models.
            # It doesn't do any harm to keep, though, so keep it for the time being.

            if ( exists $data_ref->{'wanted_gene'}->{$gene} ) { 
                # Avoid having real data overwritten by bogus '0.00' data, from split CDSes.
                if ( exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$gene}->{'rsem'} ) {
                    my $prev_rsem = $data_ref->{'data_set'}->{$label}->{'gene'}->{$gene}->{'rsem'};
                    # Keep raising the $rsem until maximized; pass over submaximal values.
                    if ( $rsem > $prev_rsem ) { 
                        $data_ref->{'data_set'}->{$label}->{'gene'}->{$gene}->{'rsem'} = $rsem;
                    }
                }
                if (! exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$gene}->{'rsem'} ) {
                    $data_ref->{'data_set'}->{$label}->{'gene'}->{$gene}->{'rsem'} = $rsem;
                }

                # Define the pseudozero for a given data set as being the smallest observed non-zero $rsem within that given data set.
                # Different data sets will probably end up having different pseudozeros.
                # Unlike the previous arbitary pseudozero of '0.01', this equals the smallest observable read count per gene.
                # Note that this can and will be overriden if $pseudozero is set to >0, e.g., to 0.01.

                # I.e., if we have not explicitly set an arbitrary non-zero number as the one-size-fits-all pseudozero for all data sets.
                # For each data set, this will drive an set-specific empirical pseudozero to be the smallest observed non-zero value.
                if ( $pseudozero == 0 ) { 
                    if ( exists $data_ref->{'data_set'}->{$label}->{'pseudocount' } ) { 
                        if ( ( $rsem > 0 ) and ( $rsem < $data_ref->{'data_set'}->{$label}->{'pseudocount' } ) ) { 
                            $data_ref->{'data_set'}->{$label}->{'pseudocount' } = $rsem;
                        }
                    }
                    if ( (! exists $data_ref->{'data_set'}->{$label}->{'pseudocount' } ) and ( $rsem > 0 ) ) { 
                        $data_ref->{'data_set'}->{$label}->{'pseudocount' } = $rsem;
                    }
                }

                # I.e., if we *have* explicitly set an arbitrary non-zero number as the one-size-fits-all pseudozero for all data sets:
                if ( $pseudozero > 0 ) { 
                    $data_ref->{'data_set'}->{$label}->{'pseudocount' } = $pseudozero;
                }
            }
        }
    }
    close $DATA_FILE or die "Can't close filehandle to data file $reports_w_gene_data[$i] (labeled $labels_of_reports[$i]): $!";
    if (! exists $data_ref->{'data_set'}->{$labels_of_reports[$i]}->{'pseudocount'} ) {
        die "Failure to record pseudocount for data file $reports_w_gene_data[$i] (labeled $labels_of_reports[$i])\n";
    }
}

if (@wanted_ratios) { 
    foreach my $ratio_name (@wanted_ratios) { 
        if ( exists $data_ref->{'data_set'}->{$ratio_name} ) { 
            die "Redundant ratio name $ratio_name\n";
        }
        if ( $ratio_name !~ /\A [^\s\/]+ \/ [^\s\/]+ \z/xms ) { 
            die "Can't parse ratio: $ratio_name\n";
        }
        if ( $ratio_name =~ /\A ([^\s\/]+) \/ ([^\s\/]+) \z/xms ) {
            my $numerator_label   = $1;
            my $denominator_label = $2;
            if ( (! exists $data_ref->{'data_set'}->{$numerator_label} ) or (! exists $data_ref->{'data_set'}->{$denominator_label} ) ) { 
                die "For ratio $ratio_name, can't find both $numerator_label and $denominator_label in labels of data sets\n";
            }
            foreach my $gene_of_int (@genes_of_interest) {
                my $numerator   = $data_ref->{'data_set'}->{$numerator_label}->{'gene'}->{$gene_of_int}->{'rsem'};
                my $denominator = $data_ref->{'data_set'}->{$denominator_label}->{'gene'}->{$gene_of_int}->{'rsem'};

                # Note that we need non-zero pseudocounts for both $numerator and $denominator, because the log10 will choke on zero for either!
                if ( $numerator == 0 ) {
                    if (! exists $data_ref->{'data_set'}->{$numerator_label}->{'pseudocount'} ) {
                        die "Numerator data set $numerator_label has no recorded pseudocount\n";
                    }
                    $numerator = $data_ref->{'data_set'}->{$numerator_label}->{'pseudocount'};
                }
                if ( $denominator == 0 ) { 
                    if (! exists $data_ref->{'data_set'}->{$denominator_label}->{'pseudocount'} ) { 
                        die "Denominator data set $denominator_label has no recorded pseudocount\n";
                    }
                    $denominator = $data_ref->{'data_set'}->{$denominator_label}->{'pseudocount'};
                }

                my $ratio      = ( $numerator  / $denominator );
                my $log10ratio = ( log($ratio) / log(10)   );
                $ratio         = sprintf("%.4f", $ratio);
                $log10ratio    = sprintf("%.4f", $log10ratio);
                $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'ratio'}      = $ratio;
                $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'log10ratio'} = $log10ratio;
            }
        }
    }
}

# Print header line:
print "Gene";

# Use grep to get rid of empty values:
@labels_of_reports = grep { $_ =~ /\S/ } @labels_of_reports;
@wanted_ratios     = grep { $_ =~ /\S/ } @wanted_ratios;

# Print actual data columns:
foreach my $label1 (@labels_of_reports) { 
    print "\t", $label1, ;
}

# Print non-zero data columns:
foreach my $label1 (@labels_of_reports) {
    print "\t", $label1, '_nz';
}

# Print ratios.
foreach my $ratio_name (@wanted_ratios) { 
    print "\t", $ratio_name, "\tlog10(", $ratio_name, ')', ;
}
print "\n";

foreach my $gene_of_int (@genes_of_interest) { 
    print "$gene_of_int";

    # Print actual data values:
    foreach my $label2 (@labels_of_reports) {
        print "\t";
        if ( exists $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rsem'} ) { 
            print $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rsem'};
        }
        else {
            print "0";  # Numerically identical to 0.00, but visibly a fill-in value; lets me spot holes in data.
        }
    }

    # Print non-zero data values:
    foreach my $label2 (@labels_of_reports) {
        print "\t";
        if (     ( exists $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rsem'} ) 
             and ( $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rsem'} > 0    ) ) { 
            print $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rsem'} ;
        }
        else { 
            print $data_ref->{'data_set'}->{$label2}->{'pseudocount'};
        }
    }

    # Print ratios, both actual and log10:
    foreach my $ratio_name (@wanted_ratios) {
        print "\t", $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'ratio'}, 
              "\t", $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'log10ratio'}, ;
    }
    print "\n";
}

