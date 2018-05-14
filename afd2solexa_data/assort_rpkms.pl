#!/usr/bin/env perl

# assort_rpkms.pl -- Erich Schwarz <emsch@caltech.edu>, 6/22/2012.
# Purpose: given WS220-style final.rpkm files from ERANGE 3.1.0, stitch them into a single report with defined ratios and no irrelevant genes.

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

my $wbgene         = q{};
my $rpkm           = q{};
my $mult_frac      = q{};
my $rpkm_non_multi = q{};
my $i              = 0;

GetOptions ( 'data=s{,}'      => \@reports_w_gene_data,
             'labels=s{,}'    => \@labels_of_reports,
             'important:s{,}' => \@reports_w_important_genes,
             'nozeros'        => \$nozeros,
             'ratios:s{,}'    => \@wanted_ratios, 
             'pseudozero=f'   => \$pseudozero,
             'help'           => \$help,   );

my $data_count  = @reports_w_gene_data;
my $label_count = @labels_of_reports;

if ( $help or (! @reports_w_gene_data) or (! @labels_of_reports) or ( $data_count != $label_count ) ) { 
    die "Format: assort_rpkms.pl\n",
        "    --data|-d         <N files, with RPKM data for genes>\n",
        "    --labels|-l       <N labels, to put atop columns of data from each listed file; if listed, must equal N; if not listed, defaults to data file names>\n",
        "    --important|-i    <labels of files whose genes we actually want to report; defaults to all N data files>\n",
        "    --nozeros|-n      <however the files listing reported genes were chosen, only list those with a nonzero rpkm in that file set>\n",
        "    --ratios|-r       <ratios that we want (if any); list as [label1]/[label2]>\n",
        "    --pseudozero|-p   [set an arbitrary, universal, fixed pseudozero such as 0.01; the default is to the smallest observed nonzero for each data set]\n",
        "    --help|-h         [print this message]\n",
        ;
}

if (! @reports_w_important_genes) { 
    @reports_w_important_genes = @reports_w_gene_data;
}

# Convert from 1-to-N to 0-to-(N-1) numbering.
$data_count--;

# This allows important files to be identified by their labels (a lot less clunky than ID by full filenames!).
foreach my $i (0..$data_count) {
    my $label        = $labels_of_reports[$i];
    my $labeled_file = $reports_w_gene_data[$i];
    $data_ref->{'file_abbrev'}->{$label} = $labeled_file;
}

foreach my $imp_gene_file_label (@reports_w_important_genes) {
    if (! exists $data_ref->{'file_abbrev'}->{$imp_gene_file_label}) { 
        die "Can't identify which important file is being referred to by important label $imp_gene_file_label\n";
    }
    my $imp_gene_file = $data_ref->{'file_abbrev'}->{$imp_gene_file_label}; 
    open my $IMP_GENES, '<', $imp_gene_file or die "Can't open file $imp_gene_file with wanted gene names: $!";
    while (my $input = <$IMP_GENES>) { 
        chomp $input;
        if ( $input =~ /\A (WBGene\d+\S*) \s+ \S+ \s+ (\S+) \s+ /xms ) {
            $wbgene = $1;
            $rpkm   = $2;
            if ( (! $nozeros) or ( $rpkm > 0 ) ) { 
                $data_ref->{'wanted_gene'}->{$wbgene} = 1;
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
        if ( $input =~ / \A (WBGene\d+\S*) \t [^\t]* \t (\d+\.\d+) \t (\d.\d{2}) \s* \z /xms) { 
            $wbgene    = $1;
            $rpkm      = $2;
            $mult_frac = $3;

            # Sanity-check numbers.
            if ( (! ( looks_like_number $rpkm ) ) or (! ( looks_like_number $mult_frac ) ) ) {
                die "In $input: RPKM $rpkm or multi/all ratio $mult_frac (or both) does not look numerical.\n";
            }
            if ( ( $mult_frac < 0 ) or ( $mult_frac > 1 ) ) { 
                die "In $input: insane multi/all ratio $mult_frac.\n";
            }

            # Get non-multi RPKM:
            $rpkm_non_multi = ($rpkm * ( 1.00 - $mult_frac ));
            $rpkm_non_multi = sprintf("%.2f", $rpkm_non_multi);

            if ( exists $data_ref->{'wanted_gene'}->{$wbgene} ) { 
                # Avoid having real data overwritten by bogus '0.00' data, from split CDSes.
                if ( exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm'} ) {
                    my $prev_rpkm = $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm'};
                    # Keep raising the $rpkm until maximized; pass over submaximal values.
                    if ( $rpkm > $prev_rpkm ) { 
                        $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm'} = $rpkm;
                    }
                }
                if (! exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm'} ) {
                    $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm'} = $rpkm;
                }

                # Do likewise with non-multiread RPKMs:
                if ( exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm_non_multi'} ) {
                    my $prev_rpkm_non_multi = $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm_non_multi'};
                    # Excelsior!
                    if ( $rpkm_non_multi > $prev_rpkm_non_multi ) {
                        $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm_non_multi'} = $rpkm_non_multi;
                    }
                }
                if (! exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm_non_multi'} ) {
                    $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene}->{'rpkm_non_multi'} = $rpkm_non_multi;
                }

                # Define the pseudozero for a given data set as being the smallest observed non-zero $rpkm within that given data set.
                # Different data sets will probably end up having different pseudozeros.
                # Unlike the previous arbitary pseudozero of '0.01', this equals the smallest observable read count per gene.
                # Note that this can and will be overriden if $pseudozero is set to >0, e.g., to 0.01.

                if ( $pseudozero == 0 ) { 
                    if ( exists $data_ref->{'data_set'}->{$label}->{'pseudocount' } ) { 
                        if ( ( $rpkm > 0 ) and ( $rpkm < $data_ref->{'data_set'}->{$label}->{'pseudocount' } ) ) { 
                            $data_ref->{'data_set'}->{$label}->{'pseudocount' } = $rpkm;
                        }
                    }
                    if ( (! exists $data_ref->{'data_set'}->{$label}->{'pseudocount' } ) and ( $rpkm > 0 ) ) { 
                        $data_ref->{'data_set'}->{$label}->{'pseudocount' } = $rpkm;
                    }
                }
                if ( $pseudozero > 0 ) { 
                    $data_ref->{'data_set'}->{$label}->{'pseudocount' } = $pseudozero;
                }
            }
        }
        # Enforce detection of unparsed lines.
        else { 
            if ( $input ne "#gene\tlen_kb\tRPKM\tmulti\/all" ) {
                # Don't die, because there may be sporadic bad lines in mostly good data.
                warn "From data file $reports_w_gene_data[$i] and expression type $label, can't parse input line:\n",
                     "$input\n",
                  ;
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
                my $numerator   = $data_ref->{'data_set'}->{$numerator_label}->{'gene'}->{$gene_of_int}->{'rpkm'};
                my $denominator = $data_ref->{'data_set'}->{$denominator_label}->{'gene'}->{$gene_of_int}->{'rpkm'};

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
                $ratio         = sprintf("%.2f", $ratio);
                $log10ratio    = sprintf("%.2f", $log10ratio);
                $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'ratio'}      = $ratio;
                $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'log10ratio'} = $log10ratio;
            }
        }
    }
}

# Print header line:
print "Gene";

# Print actual data columns:
foreach my $label1 (@labels_of_reports) { 
    print "\t", $label1, ;
}

# Print non-zero data columns:
foreach my $label1 (@labels_of_reports) {
    print "\t", $label1, '_nz';
}

# Print non-multiread data columns:
foreach my $label1 (@labels_of_reports) {
    print "\t", $label1, '_nm';
}

# At several points below, add grep to get rid of empty values.

# Print ratios.
@wanted_ratios = grep { $_ =~ /\S/ } @wanted_ratios;
foreach my $ratio_name (@wanted_ratios) { 
    print "\t", $ratio_name, "\tlog10(", $ratio_name, ')', ;
}
print "\n";

foreach my $gene_of_int (@genes_of_interest) { 
    print "$gene_of_int";

    # Print actual data values:
    @labels_of_reports = grep { $_ =~ /\S/ } @labels_of_reports;
    foreach my $label2 (@labels_of_reports) {
        print "\t";
        if ( exists $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm'} ) { 
            print $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm'};
        }
        else {
            print "0";  # Numerically identical to 0.00, but visibly a fill-in value; lets me spot holes in data.
        }
    }

    # Print non-zero data values:
    foreach my $label2 (@labels_of_reports) {
        print "\t";
        if (     ( exists $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm'} ) 
             and ( $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm'} > 0    ) ) { 
            print $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm'} ;
        }
        else { 
            print $data_ref->{'data_set'}->{$label2}->{'pseudocount'};
        }
    }
    # Print non-multiread data values:
    foreach my $label2 (@labels_of_reports) {
        print "\t"; 
        if ( exists $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm_non_multi'} ) {
            print $data_ref->{'data_set'}->{$label2}->{'gene'}->{$gene_of_int}->{'rpkm_non_multi'};
        }
        else {
            print "0";  # Again, this flags holes in data.
        }
    }

    # Print ratios, both actual and log10:
    foreach my $ratio_name (@wanted_ratios) {
        print "\t", $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'ratio'}, 
              "\t", $data_ref->{'data_set'}->{$ratio_name}->{'gene'}->{$gene_of_int}->{'log10ratio'}, ;
    }
    print "\n";
}

