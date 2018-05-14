#!/usr/bin/env perl

# read_gffs4ncbi_tbl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/6/2008.
# Purpose: map a GFF to one or more *.tbl files, for NCBI's tbl2asn.

use strict;
use warnings;

my $seq2coords_ref;
my $gene_coords_ref;
my $seq2annots_ref;

my $add_REFERENCE_line = 0;

# This is, apparently, required.  Set as coded variable, for easy changes:
my $ncbi_prefix = 'lcl|';

# Used by tsv2tbl() to extract unwanted parts of final .tbl format:
my $deletable_ref;

# Annotated during GFF reading, then used by tsv2tbl()
#     to add 3-nt stop codons to 3'-most CDSes:
my $add_stop2cds_ref;


# Read data from GFF (e.g., a Twinscan annotation of a contig):

while (my $input = <> ) { 
    chomp $input;

    # Record the range of a parent DNA sequence (*.fsa):
    if ( $input =~ / \A 
                     (\S+) 
                     \s+ assembly \s+ contig \s+ 
                     (\d+) \s+ (\d+) 
                   /xms ) { 
        my $seq1      = $1;
        my $start_nt1 = $2;
        my $end_nt1   = $3;
        $seq2coords_ref->{$seq1}->{'range'} = [$start_nt1, $end_nt1];
    }

    # Record salient gene features:
    if ( $input =~ / \A 
                     (\S+) 
                     \s+ \S+ \s+ 
                     ((?:start|stop)_codon|CDS) \s+ 
                     (\d+) \s+ (\d+) 
                     \s+ \S+ \s+ 
                     (\+|-) 
                     \s+ \S+ \s+ gene_id \s+ \" 
                     ([^\"]+) 
                     \"   
                   /xms ) { 
        my $seq2         = $1;
        my $feature_type = $2;
        my $start_nt2    = $3;
        my $end_nt2      = $4;
        my $orientation1 = $5;
        my $gene1        = $6; 

        # If a gene's termini are known, record that fact:
        if ( ( $feature_type eq 'start_codon' ) 
             or ( $feature_type eq 'stop_codon' ) ) { 
            $gene_coords_ref->{$seq2}->{$gene1}->{$feature_type} = 1;
        }

        # For later extension of 3'-most CDSes to include stop codons:
        if ( $feature_type eq 'stop_codon' ) { 
            if ( $orientation1 eq '+') { 
                my $end_cds_nt1 = ($start_nt2 - 1);
                $add_stop2cds_ref->{$seq2}->{$gene1}->{$end_cds_nt1} = $end_nt2;
            }
            if ( $orientation1 eq '-') { 
                my $end_cds_nt2 = ($end_nt2 + 1);
                $add_stop2cds_ref->{$seq2}->{$gene1}->{$end_cds_nt2} = $start_nt2;
            }
        }

        # Keep track of known gene extent, and of gene orientation:
        min_max_gene($seq2, $gene1, $start_nt2, $end_nt2);
        orient_gene($seq2, $gene1, $orientation1);

        # Unless this $feature_type's one that we want to censor from the .tbl output, 
        # encode it as .tsv line, to avoid Perl's limitations on arrayrefs as hash keys.
        # (Yes, I do know about Tie::RefHash.  I didn't see how to use it here, though.)

        my %unwanted_feature = ( 'start_codon' => 1,
                                 'stop_codon'  => 1, );

        if (! $unwanted_feature{$feature_type} ) { 
            my $info_key = $start_nt2
                           . "\t"
                           . $end_nt2
                           . "\t"
                           . $orientation1
                           . "\t"
                           . $feature_type
                           . "\t"
                           . 'gene'
                           . "\t"
                           . $gene1
                           ;

            $seq2coords_ref->{$seq2}->{$gene1}->{'coords'}->{$info_key} = 1;
        }
    }
}

# Convert GFF-derived data into NCBI *.tbl format.

foreach my $seq3 ( sort keys %{ $seq2coords_ref } ) { 
    # Make a separate *.tbl.part file for each sequence (*.fsa file).

    # This is '.tbl.part' because, to fully work, it needs editing that I find
    #    is much more easily done in a second round with a different Perl script,
    #    extracting needed protein information from the relevant '.pep' file.

    my $tbl_file = $seq3 . '.tbl.part';

    open my $TBL, '>', "$tbl_file" 
        or die "Cannot open *.tbl file $tbl_file: $!";

    my @collected_seq2annots = ();

    # Collect all the .tsv-format lines for a given DNA sequence:

    foreach my $gen_w_lims (sort keys %{ $gene_coords_ref->{$seq3} } ) { 
        my $lim_coords_txt = $gene_coords_ref->{$seq3}->{$gen_w_lims}->{'min'}
                             . "\t"
                             . $gene_coords_ref->{$seq3}->{$gen_w_lims}->{'max'}
                             . "\t"
                             . $gene_coords_ref->{$seq3}->{$gen_w_lims}->{'orientation'}
                             . "\t"
                             . 'gene'
                             . "\t"
                             . 'gene'
                             . "\t"
                             . $gen_w_lims
                             ;

         my @gen_features_txts = sort keys %{ $seq2coords_ref->{$seq3}->{$gen_w_lims}->{'coords'} };
         push @collected_seq2annots, $lim_coords_txt;
         push @collected_seq2annots, @gen_features_txts;
    }

    # Do an NCBI-.tbl-style sort of .tsv lines:
    my $presorted_seq2annots_ref = ncbi_jumblesort(\@collected_seq2annots, $seq3);
    my @presorted_seq2annots = @{ $presorted_seq2annots_ref };

    # Format the .tsv lines into NCBI .tbl line style.
    my @sorted_seq2annots = map { tsv2tbl($_, $seq3); } @presorted_seq2annots;

    # Finally, start printing.  
    # First, header for the entire sequence's *.tbl file.
    # Note the $ncbi_prefix (e.g., 'lcl|') put in front of the sequence ID.
    print {$TBL} ">Feature $ncbi_prefix", "$seq3 $tbl_file\n";

    # If a REFERENCE line is needed, this code is activatable if 
    #    the hand-coded $add_REFERENCE_line == '0' is reset to '1':
    if ( $add_REFERENCE_line ) { 
        my $seq_range = $seq2coords_ref->{$seq3}->{'range'};
        print {$TBL} $seq3, "\t", $seq_range->[0], "\t", $seq_range->[1], "\tREFERENCE\n";
    }

    # Then, a well-ordered *.tbl-format report (one hopes!):
    foreach my $annots_txt (@sorted_seq2annots) { 
        print {$TBL} "$annots_txt\n";
    }
}

# Keep a running tally of the lowest and highest nt for a given gene:
sub min_max_gene {
    my ( $seq_id, $gene_id, $start_nt, $end_nt ) = @_;

    if (    ( $seq_id   !~ / \A \S+ \z /xms ) 
         or ( $gene_id  !~ / \A \S+ \z /xms )
         or ( $start_nt !~ / \A \d+ \z /xms )
         or ( $end_nt   !~ / \A \d+ \z /xms ) ) { 
        die "Cannot parse input arguments: $seq_id, $gene_id, $start_nt, $end_nt\n";
    }

    if ( exists $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} ) { 
        if ( $start_nt <= $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} ) { 
            $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} = $start_nt;
        }
    }

    if (! exists $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} ) { 
        $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} = $start_nt;
    }

    if ( exists $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} ) { 
        if ( $end_nt >= $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} ) { 
            $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} = $end_nt;
        }
    }

    if (! exists $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} ) { 
        $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} = $end_nt;
    }
}

# Keep track of a gene's orientation; enforce consistency.
sub orient_gene { 
    my ($seq_id, $gene_id, $orientation) = @_;

    if (    ( $seq_id      !~ / \A \S+ \z /xms     )
         or ( $gene_id     !~ / \A \S+ \z /xms     )
         or ( $orientation !~ / \A (?:\+|-) \z /xms ) ) { 
        die "Cannot parse input arguments: $seq_id, $gene_id, $orientation\n";
    }

    if ( exists $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} ) { 
        if ( $orientation ne $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} ) { 
            die "Seq. $seq_id / gene $gene_id: apparently inconsistent orientation!\n";
        }
    }

    if (! exists $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} ) { 
         $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} = $orientation;
    }
}

# Convolute order of .tsv lines into NCBI's .tbl sorting style:
sub ncbi_jumblesort { 
    my ($input_array_ref, $seq_id) = @_;
    my $subarrays_ref;

    # Because Perl hashes choke on ref keys, @input_array klunkily stores .tsv lines:
    my @input_array     = @{ $input_array_ref };

    my %seen_gene       = ();
    my @genes_along_DNA = ();
    my @naively_sorted  = ();
    my @sorted_array    = ();

    # First, sort lines like a non-NCBI sane person, but keep as arrayrefs:
    @naively_sorted  = # (2) Exploit easy sorting on arrayrefs:
                       sort { $a->[0] <=> $b->[0] }
                       sort { $b->[1] <=> $a->[1] }

                       # (1) Convert from .tsv lines to arrayrefs:
                       map { my @values1 = split /\t/, $_;
                             my $values_ref1 = \@values1;
                             $values_ref1;
                       }
                       @input_array;

    # Get a clean list of genes ordered along the DNA,
    #     while simultaneously making subarrays for each gene.
    #     N.B. No provision made for genes-within-genes!

    foreach my $naive_lineref (@naively_sorted) { 
        # If this is a gene-spanning header line:
        if (    ( $naive_lineref->[-3] eq 'gene'  ) 
            and ( $naive_lineref->[-2] eq 'gene'  )
            and ( $naive_lineref->[-1] =~ /\S/xms ) ) { 

            # Record the gene along the DNA:
            my $naive_gene1 = $naive_lineref->[-1];
            if (! exists $seen_gene{$naive_gene1} ) { 
                $seen_gene{$naive_gene1} = 1;
                push @genes_along_DNA, $naive_gene1;
            }

            # Convert gene's span line from arrayref to .tsv line, 
            #    and then put it at start of gene's sorted subarray:
            my @gene_span_array = @{ $naive_lineref };
            my $gene_span_tsv   = join "\t", @gene_span_array;
            push @{ $subarrays_ref->{'final'}->{$naive_gene1} }, $gene_span_tsv;
        }

        # Otherwise, just store CDS (etc.) lines in pre-sorted subarray:
        if (    ( $naive_lineref->[-3] ne 'gene'  )
            and ( $naive_lineref->[-2] eq 'gene'  )
            and ( $naive_lineref->[-1] =~ /\S/xms ) ) { 
            my $naive_gene2 = $naive_lineref->[-1];
            push @{ $subarrays_ref->{'orig'}->{$naive_gene2} }, $naive_lineref;
        }
    }

    # Then, sort each gene's little subarrays, with the sort polarity
    #     depending on the gene's orientation!  Whee!
    #     Oh, and also, reformat back from arrayrefs to .tsv lines.
    foreach my $naive_gene3 (@genes_along_DNA) { 
        my @final_list = ();
        if ( $gene_coords_ref->{$seq_id}->{$naive_gene3}->{'orientation'} eq '+' ) { 
            @final_list = map { my @values2    = @{ $_ };
                                my $value_txt2 = join "\t", @values2;
                                $value_txt2;
                          }
                          sort { $a->[0] <=> $b->[0] }
                          @{ $subarrays_ref->{'orig'}->{$naive_gene3} };
        }
        if ( $gene_coords_ref->{$seq_id}->{$naive_gene3}->{'orientation'} eq '-' ) {
            @final_list = map { my @values3    = @{ $_ };
                                my $value_txt3 = join "\t", @values3;
                                $value_txt3;
                          }
                          sort { $b->[0] <=> $a->[0] }
                          @{ $subarrays_ref->{'orig'}->{$naive_gene3} };
        }
        # Remember: first sorted line is gene-span header! so push rest:
        push @{ $subarrays_ref->{'final'}->{$naive_gene3} }, @final_list;
    }

    # Finally, unload these subarrays into a usable single array...
    foreach my $naive_gene4 (@genes_along_DNA) { 
        my @subsorted_array = @{ $subarrays_ref->{'final'}->{$naive_gene4} };
        push @sorted_array, @subsorted_array;
    }
    return \@sorted_array;
}

# Massage .tsv lines into NCBI's crazy moon language:
sub tsv2tbl { 
    my ($input, $seq_id) = @_;

    # Reject inputs that don't fit the expected template:
    if ( $input !~ /\A \d+ \t \d+ \t (?:[+]|[-]) \t [^\t]+ \t [^\t]+ \t [^\t]+ \z /xms ) {
        die "Malformatted input: $input\n";
    }

    # If in '-' orientation, switch around nt values to descending, and splice out '-' column:
    if ( $input =~ /\A (\d+) \t (\d+) \t [-] \t ([^\t]+ \t [^\t]+ \t [^\t]+) \z /xms ) { 
        my $nt1A     = $1;
        my $nt2A     = $2;
        my $txtA     = $3;
        if ( $nt1A > $nt2A ) { 
            die "Apparently jumbled nt values ($nt1A and $nt2A) in: $input\n";
        }
        $input = $nt2A . "\t" . $nt1A . "\t" . $txtA;
    }

    # If in '+' orientation, sanity-check the nt values and splice out the '+' column:
    if ($input =~ /\A (\d+) \t (\d+) \t [+] \t ([^\t]+ \t [^\t]+ \t [^\t]+) \z /xms ) {
        my $nt1B     = $1;
        my $nt2B     = $2;
        my $txtB     = $3;
        if ( $nt1B > $nt2B ) {
            die "Apparently jumbled nt values ($nt1B and $nt2B) in: $input\n";
        }
        $input  = $nt1B . "\t" . $nt2B . "\t" . $txtB;
    }

    # Do two things to gene ends.
    # First, mark *incomplete* gene termini, running off the contig's end.
    # At the same time, add 3 nt to the 3'-most end of a CDS if a stop codon exists.
    if ($input =~ /\A (\d+) \t (\d+) \t ([^\t]+ \t [^\t]+ \t ([^\t]+)) \z /xms ) { 
        my $nt1C     = $1;
        my $nt2C     = $2;
        my $txtC     = $3;
        my $gene_id  = $4;

        # Need to have marked outputs, while keeping clean-number inputs unchanged:
        my $nt1D     = $nt1C;
        my $nt2D     = $nt2C;

        # Reject genes with poorly defined terminal coordinates:
        if (    (! $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'}  ) 
             or (! $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'}  )
             or (   $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} 
                 == $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} ) ) { 
            die "$gene_id apparently has no defined min./max.!\n";
        }

        # Mark *incomplete* gene ends with '<' or '>':
        if ( ( $nt1C == $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} ) 
             and ( ragged_low_end($seq_id,$gene_id) ) ) { 
            $nt1D = '<' . $nt1C;
        }
        if ( ( $nt1C == $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} ) 
             and ( ragged_high_end($seq_id,$gene_id) ) ) {
            $nt1D = '>' . $nt1C;
        }
        if ( ( $nt2C == $gene_coords_ref->{$seq_id}->{$gene_id}->{'min'} ) 
             and ( ragged_low_end($seq_id,$gene_id) ) ) {
            $nt2D = '<' . $nt2C;
        }
        if ( ( $nt2C == $gene_coords_ref->{$seq_id}->{$gene_id}->{'max'} ) 
             and ( ragged_high_end($seq_id,$gene_id) ) ) {
            $nt2D = '>' . $nt2C;
        }

        # Add 3 nt to 3'-ends of *complete* 3'-most CDSes:
        if ( exists $add_stop2cds_ref->{$seq_id}->{$gene_id}->{$nt1C} ) {
            $nt1D = $add_stop2cds_ref->{$seq_id}->{$gene_id}->{$nt1C};
        }
        if ( exists $add_stop2cds_ref->{$seq_id}->{$gene_id}->{$nt2C} ) {
            $nt2D = $add_stop2cds_ref->{$seq_id}->{$gene_id}->{$nt2C};
        }

        $input = $nt1D . "\t" . $nt2D . "\t" . $txtC;
    }

    # Strictly enforce correct input at this step:
    if (    ( $input !~ / \A (?: (?: [^\t]+ \t){2} [^\t]+ ) \t (?: [^\t]+ \t [^\t]+ ) \z /xms ) 
         or ( $input !~ / \A (?:<|>){0,1}\d+ \t (?:<|>){0,1}\d+ \t [^\t]+ \t [^\t]+ \t [^\t]+ \z /xms ) ) { 
        die "Can't NCBI-ize: $input\n";
    }

    # Delete all but the first qualifier of a 'gene', and all but the first 'CDS':
    if ( $input =~ / \A ( (?:<|>){0,1}\d+ \t (?:<|>){0,1}\d+ ) \t ( [^\t]+ ) \t ( [^\t]+ \t ([^\t]+) ) \z /xms ) { 
         my $coordinates1    = $1;
         my $feature1        = $2;
         my $qualifier_line1 = $3;
         my $gene_id2        = $4;
         $input              = $coordinates1;
         if (! exists $deletable_ref->{$seq_id}->{$gene_id2}->{$feature1} ) { 
             $deletable_ref->{$seq_id}->{$gene_id2}->{$feature1} = 1;
             $input .= ("\t" . $feature1);
         }
         if (! exists $deletable_ref->{$seq_id}->{$gene_id2}->{$qualifier_line1} ) { 
             $deletable_ref->{$seq_id}->{$gene_id2}->{$qualifier_line1} = 1;
             # Prefix $ncbi_prefix (e.g., 'lcl|') to any gene ID:
             $qualifier_line1 =~ s/gene\t([^\t]+)\z/gene\t$ncbi_prefix$1/;
             $input .= ("\t" . $qualifier_line1);
         }
    }

    # NCBI-ize input, if it still has a qualifier:
    if ( $input =~ / \A ( (?: [^\t]+ \t){2} [^\t]+ ) \t ( [^\t]+ \t [^\t]+ ) \z /xms ) {
        my $main_line      = $1;
        my $qualifier_line = $2;
        $input = $main_line . "\n" . "\t\t\t" . $qualifier_line;
    }

    return $input;

    # For some written documentation of NCBI's .tsv format, see:
    #     http://www.ncbi.nlm.nih.gov/Sequin/table.html#fig2
    #     
    # For an instance which makes it somewhat clearer, see:
    #     http://www.ncbi.nlm.nih.gov/Sequin/QuickGuide/sequin.htm#FeatureTableFormat
}

sub ragged_low_end { 
    my ($seq_id, $gene_id) = @_;
    my $value = 0;
    if (     ( $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} eq '+' ) 
         and (! $gene_coords_ref->{$seq_id}->{$gene_id}->{'start_codon'}       ) ) { 
        $value = 1;
    }
    if (     ( $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} eq '-'  )
         and (! $gene_coords_ref->{$seq_id}->{$gene_id}->{'stop_codon'}         ) ) { 
        $value = 1;
    }
    return $value;
}

sub ragged_high_end {
    my ($seq_id, $gene_id) = @_;
    my $value = 0;
    if (     ( $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} eq '+' )  
         and (! $gene_coords_ref->{$seq_id}->{$gene_id}->{'stop_codon'}       ) ) {
        $value = 1;
    }
    if (     ( $gene_coords_ref->{$seq_id}->{$gene_id}->{'orientation'} eq '-'  ) 
         and (! $gene_coords_ref->{$seq_id}->{$gene_id}->{'start_codon'}         ) ) {
        $value = 1;
    }
    return $value;
}

