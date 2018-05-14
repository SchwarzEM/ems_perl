#!/usr/bin/env perl

# mean_sd_rpkms.pl -- Erich Schwarz <emsch@caltech.edu>, 6/23/2012.
# Purpose: given N WS220-style final.rpkm files from ERANGE 3.1.0 (and, hopefully, 3.2.0), get a mean and SD from them.

use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;

my $format                    = 'plain';
my $label                     = q{};
my @reports_w_gene_data       = ();
my @reports_w_important_genes = ();
my $pseudozero                = 0;
my $nozeros;
my $help;

my $data_ref;

my $wbgene = q{};
my $rpkm   = q{};
my $i      = 0;

GetOptions ( 'format:s'       => \$format,
             'label:s'        => \$label,
             'data=s{,}'      => \@reports_w_gene_data,
             'important:s{,}' => \@reports_w_important_genes,
             'pseudozero=f'   => \$pseudozero,
             'nozeros'        => \$nozeros,
              
             'help'           => \$help,   );

my $data_count  = @reports_w_gene_data;

if (    $help 
     or (! @reports_w_gene_data                              ) 
     or ( ( $format ne 'plain' ) and ( $format ne 'erange' ) ) 
   ) { 
    die "Format: mean_sd_rpkms.pl\n",
        "    --format|-f       [either 'plain' or 'erange'; default is 'plain'; 'erange' is a WS-aware version of ERANGE 3.1.0 output]\n",
        "    --label|-l        [label for data being meaned/SDed; requires 1+ non-space characters; default is no label]\n",
        "    --data|-d         <N files, with RPKM data for genes>\n",
        "    --important|-i    <files whose genes we actually want to report; defaults to all N data files>\n",
        "    --pseudozero|-p   [set an arbitrary, universal, fixed pseudozero for data values, such as 0.01;\n",
        "                       the default, for each data set, is that set's smallest observed nonzero RPKM]\n",
        "    --nozeros|-n      <however the files listing reported genes were chosen, only list those with a nonzero rpkm in that file set>\n",
        "    --help|-h         [print this message]\n",
        ;
}

if (! @reports_w_important_genes) { 
    @reports_w_important_genes = @reports_w_gene_data;
}

if ( ( $label ne q{} ) and ( $label !~ /\S/ ) ) { 
    die "Label \"$label\" needs to be either blank or have some non-space characters.\n";
}

foreach my $imp_gene_file (@reports_w_important_genes) { 
    open my $IMP_GENES, '<', $imp_gene_file or die "Can't open file $imp_gene_file with wanted gene names: $!";
    while (my $input = <$IMP_GENES>) { 
        chomp $input;
        if (    ( ( $format eq 'plain'  ) and ( $input =~ /\A (WBGene\d+\S*) \s+ (\S+) \s* \z /xms      ) ) 
             or ( ( $format eq 'erange' ) and ( $input =~ /\A (WBGene\d+\S*) \s+ \S+ \s+ (\S+) \s+ /xms ) ) 
           ) {
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

# Convert from 1-to-N to 0-to-(N-1) numbering.
$data_count--;

foreach my $i (0..$data_count) { 
    # The 'label' may be a really klunky file name, but this lets me leave the rest of the code unchanged:
    my $label = $reports_w_gene_data[$i];

    open my $DATA_FILE, '<', $reports_w_gene_data[$i] or die "Can't open data file $reports_w_gene_data[$i]: $!";
    while (my $input = <$DATA_FILE>) { 
        chomp $input;
        if (    ( ( $format eq 'plain'  ) and ( $input =~ /\A (WBGene\d+\S*) \s+ (\S+) \s* \z /xms      ) )
             or ( ( $format eq 'erange' ) and ( $input =~ /\A (WBGene\d+\S*) \s+ \S+ \s+ (\S+) \s+ /xms ) )
           ) {
            $wbgene = $1;
            $rpkm   = $2;
            if ( exists $data_ref->{'wanted_gene'}->{$wbgene} ) { 
                $data_ref->{'data_set'}->{$label}->{'gene'}->{$wbgene} = $rpkm;

                # Define the pseudozero for a given data set as being the smallest observed non-zero $rpkm within that given data set.
                # Different data sets will probably end up having different pseudozeros.
                # Unlike the previous arbitary pseudozero of '0.01', this equals the smallest observable read count per gene.
                # Note that this can and will be overriden if $pseudozero is set to >0, e.g., to 0.01.

                if ( $pseudozero == 0 ) { 
                    if ( exists $data_ref->{'data_set'}->{$label}->{'pseudocount'} ) {
                        if ( ( $rpkm > 0 ) and ( $rpkm < $data_ref->{'data_set'}->{$label}->{'pseudocount'} ) ) {
                            $data_ref->{'data_set'}->{$label}->{'pseudocount'} = $rpkm;
                        }
                    }
                    if ( (! exists $data_ref->{'data_set'}->{$label}->{'pseudocount'} ) and ( $rpkm > 0 ) ) {
                        $data_ref->{'data_set'}->{$label}->{'pseudocount'} = $rpkm;
                    }
                }
                if ( $pseudozero > 0 ) { 
                    $data_ref->{'data_set'}->{$label}->{'pseudocount'} = $pseudozero;
                }
            }
        }
    }
    close $DATA_FILE or die "Can't close filehandle to data file $reports_w_gene_data[$i]: $!";
    if (! exists $data_ref->{'data_set'}->{$reports_w_gene_data[$i]}->{'pseudocount'} ) {
        die "Failure to record pseudocount for data file $reports_w_gene_data[$i]\n";
    }
}

foreach my $wanted_gene (sort keys %{ $data_ref->{'wanted_gene'} }) { 
    foreach my $i (0..$data_count) {
        # Again, this one klunky line keeps the rest of the code unchanged:
        my $label = $reports_w_gene_data[$i];

        if (! exists $data_ref->{'data_set'}->{$label}->{'gene'}->{$wanted_gene} ) { 
            die "No data for gene $wanted_gene in data set $label\n";
        }
        if ( $data_ref->{'data_set'}->{$label}->{'gene'}->{$wanted_gene} > 0 ) { 
            push @{ $data_ref->{'gene'}->{$wanted_gene}->{'vals'} }, $data_ref->{'data_set'}->{$label}->{'gene'}->{$wanted_gene};
        }
        if ( $data_ref->{'data_set'}->{$label}->{'gene'}->{$wanted_gene} == 0 ) { 
            push @{ $data_ref->{'gene'}->{$wanted_gene}->{'vals'} }, $data_ref->{'data_set'}->{$label}->{'pseudocount'};            
        }
    }
}

foreach my $wanted_gene (sort keys %{ $data_ref->{'wanted_gene'} }) { 
    my $data_count = @{ $data_ref->{'gene'}->{$wanted_gene}->{'vals'} };
    if ( $data_count < 3 ) { 
        die "Too few data points available for gene $wanted_gene\n";
    }
    my $stat1 = Statistics::Descriptive::Full->new();
    $stat1->add_data(@{ $data_ref->{'gene'}->{$wanted_gene}->{'vals'} });
    $data_ref->{'gene'}->{$wanted_gene}->{'mean'}    = $stat1->mean();
    $data_ref->{'gene'}->{$wanted_gene}->{'std_dev'} = $stat1->standard_deviation();
}

# Print header line:
print "Gene\t", $label, "Mean\t", $label, "SD\t", $label, "z-over-0\n", ;

foreach my $gene_of_int (@genes_of_interest) { 
    my $mean    = $data_ref->{'gene'}->{$gene_of_int}->{'mean'};
    $mean       = sprintf("%.2f", $mean);
    my $sd      = $data_ref->{'gene'}->{$gene_of_int}->{'std_dev'};
    $sd         = sprintf("%.2f", $sd);

    my $z_score = 'n/a';
    if ( $data_ref->{'gene'}->{$gene_of_int}->{'std_dev'} > 0 ) { 
        $z_score = ( $data_ref->{'gene'}->{$gene_of_int}->{'mean'} / $data_ref->{'gene'}->{$gene_of_int}->{'std_dev'} );
        $z_score = sprintf("%.2f", $z_score);
    }

    print "$gene_of_int\t$mean\t$sd\t$z_score\n";
}

