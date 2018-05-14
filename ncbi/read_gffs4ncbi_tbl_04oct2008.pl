#!/usr/bin/env perl

# read_gffs4ncbi_tbl_04oct2008.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/4/2008.
# Purpose: based on a poorly-documented file specification, mis-map a GFF to one or more *.tbl files, for NCBI's tbl2asn.

use strict;
use warnings;

my $seq2coords_ref;
my $gene_coords_ref;
my $seq2annots_ref;

my $add_REFERENCE_line = 0;

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
    # Make a separate *.tbl file for each sequence (*.fsa file):
    my $tbl_file = $seq3 . '.tbl';
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

    # Order and format the .tsv lines, with a longish sort-map:

    my @sorted_seq2annots = # [5] Then, map each .tsv line into wonky NCBI .tbl format.
                            # $seq3 has to be fed in to allow marking of ragged gene ends;
                            #    the gene ID is also needed, but is carried in with $_.
                            map { tsv2tbl($_, $seq3); } 

                            # [4] Convert the processed data back to .tsv text-lines:
                            map { my @values = @{ $_ };
                                  my $value_txt = join "\t", @values;
                                  $value_txt;
                            }

                            # [3] Then, exploit easy sorting on arrayrefs:
                            sort { $a->[0] <=> $b->[0] } 
                            sort { $b->[1] <=> $a->[1] }

                            # [2] First -- map back to single-scalar array references!
                            map { my @values = split /\t/, $_;
                                  my $values_ref = \@values;
                                  $values_ref; 
                            } 

                            # [1] Because Perl hashes choke on ref keys, 
                            #     this array had to klunkily store *.tsv text-lines:
                            @collected_seq2annots;

    # Finally, start printing.  First, header for the entire sequence's *.tbl file:
    print {$TBL} ">Feature $seq3 $tbl_file\n";

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

# Hairy subroutine to massage .tsv lines into NCBI's crazy moon language:
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

    # Mark *incomplete* gene termini, running off the contig's end:
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

        # In all cases below, only mark-up *incompletes*:
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
        $input = $nt1D . "\t" . $nt2D . "\t" . $txtC;
    }

    # Final step: convert $input from simple .tsv line to the weird 
    #     "three-columns / linefeed / three tabs, two columns"
    #     format that NCBI demands:

    # Enforce expected input.
    if ( $input !~ / \A (?: (?: [^\t]+ \t){2} [^\t]+ ) \t (?: [^\t]+ \t [^\t]+ ) \z /xms ) { 
        die "Can't NCBI-ize: $input\n";
    }

    # NCBI-ize input.
    if ( $input =~ / \A ( (?: [^\t]+ \t){2} [^\t]+ ) \t ( [^\t]+ \t [^\t]+ ) \z /xms ) {
        my $main_line      = $1;
        my $qualifier_line = $2;
        $input = $main_line . "\n" . "\t\t\t" . $qualifier_line;
    }

    return $input;
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

