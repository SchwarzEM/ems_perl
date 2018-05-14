#!/usr/bin/env perl

# parse_pecan_clw.pl -- Erich Schwarz <emsch@caltech.edu>, 10/18/2011.
# Purpose: given Michael Paulini's PECAN alignment in ClustalW format, extract coordinates of elegans sequences with 3+-fold identity to homologs.

use strict;
use warnings;
use Getopt::Long;

my $input_alignment      = q{};

my $input                = q{};
my $prev_line            = q{};

my $ready_for_stanza       = 1;
my $inside_stanza          = 0;
my $ready_for_stanza_end   = 0;
my $counting_seqs          = 0;

my $segment_text         = q{};
my $stable_buffer_length = 0;
my $new_buffer_length    = 0;

my $seqs_in_align    = 0;  
my $elegans_present  = 0; 
my $usable_alignment = 0;

my $elegans_segment_name = q{};
my $remembered_ele_name  = q{};
my $chr                  = q{};
my $start_nt             = 0;
my $end_nt               = 0;

my $non_elegans_segment_name = q{};

my $raw_sequence         = q{};
my $raw_cons_sites       = q{};

my $buffer_size          = 0;

my $data_ref;
my $help;

GetOptions ( "align=s"  => \$input_alignment,
             "buffer=i" => \$buffer_size, 
             "help" => \$help,     );

if ( $help or (! $input_alignment ) ) { 
    die "Format: parse_pecan_clw_19oct2011.pl",
        " --align|-a [input pecan.clw alignment]",
        " --buffer|-b [optional buffer nt size]",
        " --help|-h\n",
        ;
}

if ( ( $buffer_size < 0 ) or ( $buffer_size != int($buffer_size) ) ) { 
    die "Buffer size $buffer_size must be a non-negative integer!\n";
}
# Revert to default value if misspecified in arguments:
$buffer_size |= 0;

open my $ALIGN, '<', $input_alignment or die "Can't open input PECAN ClustalW alignment $input_alignment: $!";
while (my $newline = <$ALIGN>) { 
    # Keep a record of each previous line.
    chomp $newline;
    $prev_line = $input;
    $input     = $newline;

    if ( $input =~ /\A CLUSTAL [ ] W \( 1 \. 81 \) [ ] multiple [ ] sequence [ ] alignment \s* \z/xms ) { 
        # Enforce correct state:
        if (! $ready_for_stanza ) { 
            die "Can't parse start of stanza: $input\n";
        }

        # Invoke data-processing iff we want the alignment at all:
        if ($usable_alignment) { 
            analyze_export_and_clear_data($remembered_ele_name);
        }
        if (! $usable_alignment) { 
            delete $data_ref->{'segment'}->{$remembered_ele_name};
            delete $data_ref->{'homolog_segments'};
        }

        $ready_for_stanza     = 0;
        $inside_stanza        = 0;
        $counting_seqs        = 1;
        $seqs_in_align        = 0;
        $elegans_present      = 0;
        $usable_alignment     = 0;

        $segment_text         = q{};
        $stable_buffer_length = 0;
        $new_buffer_length    = 0;

        $elegans_segment_name = q{};
        $remembered_ele_name  = q{};

        $start_nt             = 0;
        $end_nt               = 0;  

        $raw_sequence         = q{};
        $raw_cons_sites       = q{};
    }
    elsif ( $input =~ /\A (\S+ \( (?:\-|\+) \) \/ \d+ \- \d+ \s+) [ACGTN-]{1,60} \z /xms ) { 
        # Empirically determine left buffer length for entire alignment stanza set.
        $segment_text = $1;
        $new_buffer_length = length $segment_text;
        if ($stable_buffer_length and ( $stable_buffer_length != $new_buffer_length) ) {
            die "Inconsistent buffer lengths: $stable_buffer_length vs. $new_buffer_length\n";
        }
        if (! $stable_buffer_length) {
            $stable_buffer_length = $new_buffer_length;
        }

        $inside_stanza = 1;

        if ($counting_seqs) { 
            $seqs_in_align++;
        }
        if ( $input =~ /\A Cele\S* /xms ) { 
            # Require the following rather specific elegans format, on pain of death:
            if ( $input =~ /\A (Cele\-CHROMOSOME_(I|II|III|IV|V|X) \(\+\) \/ (\d+) \- (\d+) \s+) ([ACGTN-]{1,60}) \z /xms ) { 
                $elegans_segment_name = $1;
                $chr                  = $2;
                $start_nt             = $3;
                $end_nt               = $4;
                $raw_sequence         = $5;

                # Make absolutely sure the elegans text fits the alignment's buffer length:
                $new_buffer_length = length $elegans_segment_name;
                if ( $stable_buffer_length != $new_buffer_length ) {
                    die "Mysteriously inconsistent buffer lengths: $stable_buffer_length vs. $new_buffer_length\n";
                }
                # Then trim off trailing spaces to get a simple name for the elegans aligned block:
                $elegans_segment_name =~ s/\s+\z//;

                # Enforce consistency there, too:
                if ( $remembered_ele_name and ( $remembered_ele_name ne $elegans_segment_name ) ) {
                    die "Two different elegans segments ($remembered_ele_name vs. $elegans_segment_name) in one stanza: $input\n";
                }
                if (! $remembered_ele_name) {
                    $remembered_ele_name  = $elegans_segment_name;
                }

                $elegans_present = 1;

                $data_ref->{'segment'}->{$elegans_segment_name}->{'start'}   = $start_nt;
                $data_ref->{'segment'}->{$elegans_segment_name}->{'end'}     = $end_nt;
                $data_ref->{'segment'}->{$elegans_segment_name}->{'rawseq'} .= $raw_sequence;
            }
            else { 
                die "Putative elegans line in unacceptable format: $input\n";
            }
        }
        # Record homolog rawseqs, so that I can later determine 3+-way identity directly from the data,
        #      without relying on totally flakely ClustalW starring:
        if ( $input !~ /\A Cele\S* /xms ) { 
            if ( $input =~ /\A ( \S+ \( (?:\-|\+) \) \/ \d+ \- \d+) \s+ ([ACGTN-]{1,60}) \z /xms ) {
                $non_elegans_segment_name = $1;
                $raw_sequence             = $2;
                $data_ref->{'homolog_segments'}->{$non_elegans_segment_name}->{'rawseq'} .= $raw_sequence;
            }
            else { 
                die "Can't extract non-elegans data out of: $input\n";
            }
        }
    }
    elsif ( $input =~ / \A \s{$stable_buffer_length} (.+) \z /xms ) {
        $raw_cons_sites       = $1;
        $counting_seqs        = 0;
        $ready_for_stanza_end = 1;
        if ( $elegans_present and ( $seqs_in_align >= 3 ) ) {
            $usable_alignment = 1;
            $data_ref->{'segment'}->{$elegans_segment_name}->{'rawcons'} .= $raw_cons_sites;
        }
        else {
            $usable_alignment = 0;
            delete $data_ref->{'segment'}->{$elegans_segment_name};
        }
    }
    elsif ( $input !~ /\A .+ \z /xms ) { 
        if ($inside_stanza) { 
            if ($ready_for_stanza_end) { 
                $ready_for_stanza_end = 0;
                $inside_stanza        = 0;
                $ready_for_stanza     = 1;
            }
            else { 
                die "Not ready for stanza end -- $ready_for_stanza_end\n",
                    "Current input \"$input\"\n",
                    "Past remembered line: \"$prev_line\"\n",
                    ;
            }
        }
    }
    else { 
        die "Can't parse data\n",
            "Prev. rememb. line: $prev_line\n",
            "Current input line: $input\n",
            ;
    }
}

# EOF data handling:
if ($usable_alignment) {
    analyze_export_and_clear_data($remembered_ele_name);
}

# EOF diagnostic message:
if (! $usable_alignment) {
    delete $data_ref->{'segment'}->{$remembered_ele_name};
    delete $data_ref->{'homolog_segments'};
}
close $ALIGN or die "Can't close filehandle to PECAN ClustalW alignment $input_alignment: $!";

sub analyze_export_and_clear_data {
    my $ele_name = $_[0];

    # If fed an empty 'name', end subroutine silently at once:
    if ( $ele_name !~ /\S/xms ) { 
        return;
    }

    # Get homolog list and make sure it has at least two members:
    my @homolog_list = sort keys %{ $data_ref->{'homolog_segments'} };
    my $homolog_list_text = join '; ', @homolog_list;
    my $homolog_count = q{};
    $homolog_count = @homolog_list;
    if ( $homolog_count < 2 ) {
        die "Not enough homologs for elegans block $ele_name are observed (only $homolog_count)\n";
    }

    if ( $ele_name and ( exists $data_ref->{'segment'}->{$ele_name} ) ) {
        # Ensure that initial lengths of raw sequence and alignment hit lines are consistent.
        # I don't actually use this for scoring alignment any more becuase it's highly unreliable,   
        #     but enforcing that I read it properly is one way to enforce a successful file parse.

        my $raw_elegans_length  = length $data_ref->{'segment'}->{$ele_name}->{'rawseq'};
        my $rawcons_length      = length $data_ref->{'segment'}->{$ele_name}->{'rawcons'};
        if ( $raw_elegans_length != $rawcons_length ) { 
            die "Raw seq. length, $raw_elegans_length, does not equal raw cons. length, $rawcons_length\n";
        }

        # Get Perl 0-based length of the raw sequence lines, then array-ize the sequences so that I can work through them.
        $raw_elegans_length--;
        my @raw_residues  = split //, $data_ref->{'segment'}->{$ele_name}->{'rawseq'};
        foreach my $hom_seq (@homolog_list) { 
            @{ $data_ref->{'segment'}->{$ele_name}->{'homologs'}->{$hom_seq}->{'indiv_res'} } 
                = split //, $data_ref->{'homolog_segments'}->{$hom_seq}->{'rawseq'}
        }

        my $j          = 0;
        # First, decide which residues in elegans are matches for 2+ homologs, and mark them as upper-case.
        # At this stage, we want to consider all residues including gaps, because otherwise we fail!
        foreach my $i (0..$raw_elegans_length) { 
            if ( $raw_residues[$i] !~ /\A[ACGTN-]\z/xms ) {
                die "Can't parse residue: $raw_residues[$i]\n";
            }
            if ( $raw_residues[$i] =~ /\A[ACGTN-]\z/xms ) { 
                my $input_res  = $raw_residues[$i];
                my $cons_score = 0;

                # Count the number of identities:
                foreach my $hom_seq1 (@homolog_list) {
                    if ( $data_ref->{'segment'}->{$ele_name}->{'homologs'}->{$hom_seq1}->{'indiv_res'}->[$i] eq $raw_residues[$i] ) {
                        $cons_score++;
                    }
                    else { 
                        if ( $data_ref->{'segment'}->{$ele_name}->{'homologs'}->{$hom_seq1}->{'indiv_res'}->[$i] !~ /\A[ACGTN-]\z/xms ) {
                            die "Illegal residue $data_ref->{'segment'}->{$ele_name}->{'homologs'}->{$hom_seq1}->{'indiv_res'} in homolog $hom_seq1\n";
                        }
                    }
                }

                # Only bother adding residues to the final sequence, or markers to the final consensus output, if they match elegans residues, not gaps!
                if ( $input_res =~ /\A[ACGTNacgtn]\z/xms ) {
                    $j++;
                    if ( $cons_score <= 1 ) { 
                        $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{$j} = 'non_conserved';
                    }
                    elsif ( $cons_score >= 2 ) {
                        $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{$j} = 'conserved';
                    }
                    else {
                        die "Can't parse cons. score $cons_score of $ele_name\n";
                    }
                }
            }
        }

        # Having worked through sequences to get final ones, process:
        my $correct_length
            = ( $data_ref->{'segment'}->{$ele_name}->{'end'}
                    - $data_ref->{'segment'}->{$ele_name}->{'start'}
                    + 1 );

        my @local_coords = sort { $a <=> $b } keys %{ $data_ref->{'segment'}->{$ele_name}->{'local_coords'} };

        # Before actually exporting data, but after we *have* the data (and thus can avoid creating imaginary residues),
        #    push outward from the conserved residues *if* we have set a buffer size.
        #    Note that we do not want to cause an autocatalytic spread of 'conserved' markers, so we use a 'buffered' tag that is inert to further change.
        if ($buffer_size) {
            foreach my $d (@local_coords) {
                if ( $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{$d} eq 'conserved' ) {
                    foreach my $b (1..$buffer_size) {
                        if (      ( exists $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{($d-$b)}             ) 
                              and ( $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{($d-$b)} eq 'non_conserved' ) ) { 
                            $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{($d-$b)} = 'buffered';
                        }
                        if (     ( exists $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{($d+$b)}             ) 
                             and ( $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{($d+$b)} eq 'non_conserved' ) ) { 
                            $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{($d+$b)} = 'buffered';
                        }
                    }
                }
            }
        }

        # After doing whatever buffering we need, do safety check:
        my $local_len = @local_coords;
        if ( $local_len != $correct_length ) {
            die "Local alignment length, $local_len, does not equal correct length, $correct_length\n";
        }

        my $coord_shift = ($start_nt - 1);
        my $block_start = 0;
        my $block_stop  = 0;

        foreach my $k (@local_coords) { 
             if ( $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{$k} eq 'non_conserved' ) {
                 if ($block_start) { 
                     $block_start += $coord_shift;
                     $block_stop  += $coord_shift;
                     print "$chr\t$block_start\t$block_stop\n";
                     $block_start = 0;
                     $block_stop  = 0;
                 }
             }
             elsif (    ( $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{$k} eq 'conserved' ) 
                     or ( $data_ref->{'segment'}->{$ele_name}->{'local_coords'}->{$k} eq 'buffered'  ) ) { 
                 if ($block_start) { 
                     $block_stop = $k;
                 }
                 else { 
                     $block_start = $k;
                     $block_stop  = $k;
                 }
             }
             else { 
                 die "Can't parse local residue $k of $ele_name\n";
             }
        }
        # EOF clearout:
        if ($block_start) {
            $block_start += $coord_shift;
            $block_stop  += $coord_shift;
            print "$chr\t$block_start\t$block_stop\n";
            $block_start = 0;
            $block_stop  = 0;
        }

        delete $data_ref->{'segment'}->{$ele_name};
        delete $data_ref->{'homolog_segments'};
    }
}

