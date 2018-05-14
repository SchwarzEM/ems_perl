#!/usr/bin/env perl

# merge_LCgene_class_lists.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/24/2010.
# Purpose: given six different LC gene lists, merge them into a single coherent list.  LEGACY version for earlier kludge.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;

my $gene_data_ref;
my @input_files = @ARGV;

# Hand-code acceptable inputs.  Klunky but OK for this specific job.

my %req_inputs2desc = ( 'L3_class_A.rel_rpkm.tsv' => 'class A, L3',
                        'L4_class_A.rel_rpkm.tsv' => 'class A, L4',
                        'L3_class_B.rel_rpkm.tsv' => 'class B, L3',
                        'L4_class_B.rel_rpkm.tsv' => 'class B, L4',
                        'L3_class_C.rel_rpkm.tsv' => 'class C, L3',
                        'L4_class_C.rel_rpkm.tsv' => 'class C, L4', );

if (! @input_files ) { 
    die "Format: merge_LCgene_class_lists.pl  [input files:] L(3|4)_class_(A|B|C).rel_rpkm.tsv  <STDOUT>\n";
}

foreach my $infile (@input_files) { 
    # Put these declarations inside a loop, to seal them off from 'variable carryover'.
    my $gene = q{};
    my $rpkm = 0;
    my $annot = q{};

    # Enforce acceptability.
    if (! exists $req_inputs2desc{$infile} ) { 
        die "Format: merge_LCgene_class_lists.pl  [input files:] L(3|4)_class_(A|B|C).rel_rpkm.tsv  <STDOUT>\n";
    }

    # For each data file, get genes listed, enforce consistent annotations, and gather gene RPKMs.
    open my $INFILE, '<', $infile or die "Can't open input file $infile: $!";
    while (my $input = <$INFILE>) { 
        chomp $input;
        if ( $input =~ /\A (WBGene\S+) \t \d+\.\d+ \t ([^\t]+) \t (.*) \z  /xms ) { 
            $gene  = $1;
            $rpkm  = $2;
            $annot = $3;
            if (! ( looks_like_number $rpkm ) ) {
                die "RPKM $rpkm does not look numerical.\n";
            }
            if (      ( exists $gene_data_ref->{$gene}->{'annot'}    ) 
                  and ( $annot ne $gene_data_ref->{$gene}->{'annot'} ) ) { 
                die "Inconsistent annotations for gene $gene:\n\n",
                    "$gene_data_ref->{$gene}->{'annot'}\n",
                    "$annot\n",
                    "\n",
                    ;
            }
            push @{ $gene_data_ref->{$gene}->{'RPKMs'} }, $rpkm;
            $gene_data_ref->{$gene}->{'annot'} = $annot;
            $gene_data_ref->{$gene}->{'classes'}->{ $req_inputs2desc{$infile} } = 1;
        }
        else {
            die "Can't parse input line: $input\n";
        }
    }
    close $INFILE or die "Can't close filehandle to input file $infile: $!";
}

# Get summaries of gene data:
foreach my $wbgene1 ( sort keys %{ $gene_data_ref } ) { 

    # Get the maximum observed RPKM in L3 and L4 (there should be, at most, one RPKM each):
    my $stat = Statistics::Descriptive::Sparse->new();
    $stat->add_data( @{ $gene_data_ref->{$wbgene1}->{'RPKMs'} } );
    $gene_data_ref->{$wbgene1}->{'max RPKM'} = $stat->max();
    undef $stat;

    # Get a single text string listing all class memberships; again, at most there should be two.
    $gene_data_ref->{$wbgene1}->{'class summary'}
        = join q{; }, (sort keys %{ $gene_data_ref->{$wbgene1}->{'classes'} } );

}

# For the final summary, list the genes in descending mean expression (RPKM).
foreach my $wbgene2 ( sort {     $gene_data_ref->{$b}->{'max RPKM'} 
                             <=> $gene_data_ref->{$a}->{'max RPKM'} } 
                      keys %{ $gene_data_ref } ) { 
    print $wbgene2,
          "\t",
          $gene_data_ref->{$wbgene2}->{'max RPKM'},
          "\t",
          "\"$gene_data_ref->{$wbgene2}->{'class summary'}\"",
          "\t",
          $gene_data_ref->{$wbgene2}->{'annot'},
          "\n",
          ;
}

