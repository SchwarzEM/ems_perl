#!/usr/bin/env perl

# split_fasta.pl -- Erich Schwarz <ems394@cornell.edu>, 7/22/2013.

# Purpose: given a FASTA sequence file, either split into single sequences named after themselves, or into N sub-FASTA blocks; blocks can have ~equal numbers of seqs., ~equal number of *contigs* (where we are splitting scaffolds; this use is for dealing with NCBI's defective uploading limits), or ~equal sizes.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my @input_files = ();

my $block_count   = 0;
my $groups_count  = 0;
my $contigs_count = 0;
my $split_all;

my $sizes_equal;
my $numbers_equal;

my $seq_count     = 0;
my $gr_blk_count  = 0;

my $output_name  = q{};
my $output_line  = q{};
my @output_lines = (); 
my $seq_name     = q{};

my $data_ref;

my @chosen_seqs = ();

my $help;

GetOptions ( 'input_files=s{,}' => \@input_files,
             'output:s'         => \$output_name,

             # May only select only one of these four options, or must choose only 'all':
             'blocks=i'         => \$block_count,
             'groups=i'         => \$groups_count,
             'contigs=i'        => \$contigs_count,
             'all'              => \$split_all,

             # May only select one of these two options, and only if 'blocks' is also chosen:
             'sizes'            => \$sizes_equal,
             'numbers'          => \$numbers_equal,

             'help'             => \$help, );

if (    $help 
     # Enforce logic of inputs.
     # Obviously, need input files!
     or (! @input_files )

     # Must have exactly one of these: $block_count, $groups_count, $contigs_count, or $split_all.
     or ( (! $block_count) and (! $groups_count ) and (! $contigs_count ) and (! $split_all) )
     or ( $block_count   and ( $groups_count or $contigs_count or $split_all ) ) 
     or ( $groups_count  and ( $block_count  or $contigs_count or $split_all ) )
     or ( $contigs_count and ( $block_count  or $groups_count  or $split_all ) )
     or ( $split_all     and ( $block_count  or $groups_count  or $contigs_count ) )

     # Must have $sizes_equal or $numbers_equal for $block_count, only.
     or ( $block_count and ( (! $sizes_equal ) and (! $numbers_equal) ) )
     or ( ( $groups_count or $contigs_count or $split_all ) and ( $sizes_equal or $numbers_equal ) )

     # May only have one of these:
     or ( $sizes_equal and $numbers_equal )
 
   ) {
    die "\n",
        "Format: split_fasta_sizewise.pl\n",
        "        --input_files|-i  [input file(s), or '-' if stream]\n",
        "        --output|-o       [optional: base name for output files (default is first input filename, or \"output\" if input is '-' stream)]\n",
        "\n",
        "        [Must either select only one of these four options:]\n",
        "        --blocks|-b <N>   [split into N ~equal blocks (equal either by --sizes or by --numbers)]\n",
        "        --groups|-g <N>   [split into blocks with N sequences per block ('group'); instead of controlling block count, we control top limit of membership]\n",
        "        --contigs|-c <N>  [split into blocks with N *contigs* per block; as with 'group', we control top limit of membership]\n",
        "        --all|-a          [split all individual members into self-named seqs]\n",
        "\n",
        "        [May only select one of these two options, and only if --blocks is also chosen]\n",
        "        --sizes|-s        [split the N blocks into ~equal sizes, as measured by residues (nt or aa) per block]\n",
        "        --numbers|-n      [split the N blocks into ~equal numbers of sequences per block]\n",
        "\n",
        "        --help|-h         [print this message]\n",
        "\n",
        ;
}

if (! $output_name ) { 
    if ( ( $input_files[0] ne '-') and ( $input_files[0] =~ /\A \S+ \z/xms ) ) { 
        $output_name = basename($input_files[0]);
    } 
    else { 
        $output_name = 'output';
    }
}

if ($block_count) {
    check_count_val($block_count, 'block')
}
if ($groups_count) {
    check_count_val($groups_count, 'groups')
}
if ($contigs_count) {
    check_count_val($contigs_count, 'contigs')
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
foreach my $infile (@input_files) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }

    # Record the incoming FASTA data.
    while (my $input_line = <$INPUT_FILE>) {
        chomp $input_line;   
        if ($input_line =~ /\A > ( (\S+) .*) /xms) {
            my $header = $1;
            $seq_name  = $2;
            # Enforce nonredundant names for each sequence:
            if (exists $data_ref->{'seq_name'}->{$seq_name}) { 
                die "Redundant sequence name: $seq_name\n";
            }
            $data_ref->{'seq_name'}->{$seq_name}->{'header'} = $header;

            # For a variety of different splits, keep the original order of sequence names:
            if ($groups_count or $contigs_count or $numbers_equal) { 
                push @{ $data_ref->{'original_seq_list'} }, $seq_name;
            }
        }
        elsif ( $input_line =~ /\S/xms ) {
            $data_ref->{'seq_name'}->{$seq_name}->{'seq'} .= $input_line;
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!\n";
}

# After having imported all sequence data, parse it out either into various block types, 
#     or into individual sequences uniquely named after themselves, in a new directory. 

if ($split_all) {
    foreach my $seq_id ( sort keys %{ $data_ref->{'seq_name'} } ) { 
        # Each sequence gets its own safely-named .fa (or .fa.\d+, but that's the cost of safety):
        my $seq_file = $seq_id . '.fa';
        $seq_file    = safename($seq_file);

        open my $SEQ_ID, '>', $seq_file or die "Can't open sequence output file $seq_file: $!";

        # Header:
        print $SEQ_ID '>', $data_ref->{'seq_name'}->{$seq_id}->{'header'}, "\n";

        # FASTA-formatted sequence text:
        my $curr_seq = $data_ref->{'seq_name'}->{$seq_id}->{'seq'};
        @output_lines = unpack( "a60" x ( length($curr_seq)/60 + 1), $curr_seq );
        foreach $output_line (@output_lines) {
                print $SEQ_ID "$output_line\n" if ($output_line =~ /\S/);
        }

        close $SEQ_ID or die "Can't close filehandle to sequence output file $seq_file: $!";
    }
    exit;
}

# The logic isn't very different for the various ways of allocating blocks;
#     the only thing that ends up really mattering is how sequences get mapped *to* blocks.

if ($block_count or $groups_count or $contigs_count) { 
    my $total_sizes    = 0;
    my $sum_sizes      = 0;
    my $assigned_block = 0;

    # List names:
    my @seq_names = ();

    # If block sizes are to be equal, then list names by decreasing size:
    if ($sizes_equal) { 
        @seq_names = sort {     length($data_ref->{'seq_name'}->{$b}->{'seq'})
                           <=> length($data_ref->{'seq_name'}->{$a}->{'seq'}) }
                    keys %{ $data_ref->{'seq_name'} };
    }

    # For various other sortings, recover the original order of sequences.
    #     Note: could instead sort ASCIIbetically, with 'sort keys %{ $data_ref->{'seq_name'} }':

    if ($groups_count or $contigs_count or $numbers_equal) {
        @seq_names = @{ $data_ref->{'original_seq_list'} };
    }

    # Get total sequence length.  Here, it really does matter how we are going to compute 'length'!
    foreach my $seq_name (@seq_names) { 
        # For equal sizes we need the total sequence size.
        if ($sizes_equal) {
            $total_sizes += length($data_ref->{'seq_name'}->{$seq_name}->{'seq'});
        }

        # For equal numbers or numbered groups, we will give each sequence a nominal 'size' of 1 nt!
        # Which lets me keep the count machinery for equal sizes unchanged, yet get equal-member results.
        if ($groups_count or $numbers_equal) { 
            $total_sizes = @seq_names;
        }

        # For counting out by contigs, we need to compute contigs for each scaffold, something we have not done until now.
        if ($contigs_count) {
            my @scaf_contigs   = split /[Nn]+/, $data_ref->{'seq_name'}->{$seq_name}->{'seq'};
            my $scaf_contig_no = @scaf_contigs;

            # Note that we have no really satisfactory way to deal with scaffolds which have more contigs that our theoretical limit.
            # So, we sweep them under the rug, by giving "over-contiged" scaffolds the nominal maximum allowed number of contigs.
            # We also *record* this fact, so that we have some idea what problems to expect later.
            if ( $scaf_contig_no > $contigs_count ) {
                warn "For sequence $seq_name, actually observed $scaf_contig_no contigs; nominally rounded down to $contigs_count\n";
                $scaf_contig_no = $contigs_count;
            }
            $data_ref->{'seq_name'}->{$seq_name}->{'contig_count'} = $scaf_contig_no;
            $total_sizes += $scaf_contig_no;
        }
    }

    # Assign each name to blocks from 1 to $block_count.
    # For three of our options ($contigs_count, $sizes_equal, or $numbers_equal), this depends on block-size logic. 
    # For $groups_count, we have to just count heads.

    foreach my $seq_name (@seq_names) {
        if ($sizes_equal) {
            $sum_sizes += length($data_ref->{'seq_name'}->{$seq_name}->{'seq'});
        }
        if ($numbers_equal or $groups_count) {
            # Again, to get split by numbers of sequences, just make each sequence 'length' 1 nt.
            $sum_sizes++;
        }
        if ($contigs_count) {
            # Note that we have already bowlderized 'over-contiged' counts.
            my $scaf_contig_no = $data_ref->{'seq_name'}->{$seq_name}->{'contig_count'};
            $sum_sizes += $scaf_contig_no;
        }

        # From here on, the machinery I need to get any kind split is almost exactly the same.
        # Any which way, we need to get a rough value of the block to which the sequence will be assigned.

        # For methods that split files into a predetermined number of blocks:
        if ($sizes_equal or $numbers_equal) { 
            # Getting this right requires, first, getting a rough value of the block:
            $assigned_block = ( ( $sum_sizes * $block_count ) / $total_sizes );
        }

        # For the method that allows any number of files, but enforces <= N sequences per file:
        if ($groups_count) {
            # Getting this right requires, first, getting a rough value of the block:
            # Note that we need to -1 $sum_sizes, or our *first* block will be 1 sequence short of its intended total.
            $assigned_block = ( ($sum_sizes - 1) / $groups_count);
        }

        # For the method that allows any number of files, but tries to enforce <= N *contigs* per file (while overlooking impossibilities):
        if ($contigs_count) {
            $assigned_block = ( ($sum_sizes - 1) / $contigs_count);
        }

        # After that, the logic will be almost the same.
        # Round the initial block value (yes, this is not ideal -- sprintf is supposed to be better -- but it will work):
        $assigned_block = int($assigned_block);

        # In almost all cases, what I just did rounds from 0 to $block_count-1, so add 1:
        $assigned_block++;

        # Correct for edge effects.  We should *never* accidentally get a block less than 1:
        if ( $assigned_block < 1 ) { 
            $assigned_block = 1;
        }

        # If we are splitting by N sequences or contigs per block (and not by equal blocks)
        #     then we need to keep a continuously ascending tally of what our final block number will be:
        if ( ( $groups_count or $contigs_count ) and ( $assigned_block > $gr_blk_count ) ) { 
            $gr_blk_count = $assigned_block;
        }

        # If we are splitting by equal blocks (and not by N sequences per block),
        #     then we should not have a block greater than $block_count itself, either:
        if ( (! $groups_count) and (! $contigs_count) and ( $assigned_block > $block_count ) ) { 
            $assigned_block = $block_count;
        }

        # Finally, we have a well-chosen block integer for each sequence:
        $data_ref->{'seq_name'}->{$seq_name}->{'block'} = $assigned_block ;
    } 

    # At last, for each block, print out a safely-named FASTA compendium of its members.
    # Print out sequences in order of descending size.

    # If we have been working by groups with N sequences or contigs, then put the ascended value $gr_blk_count back into $block_count:
    if ($groups_count or $contigs_count) {
        $block_count = $gr_blk_count ;
    }

    foreach my $block (1..$block_count) { 
        # Make the numbers of blocks sort well by ASCII:
        my $blk_no    = $block;
        my $DIGITS    = length($block_count);
        my $sf_format = '%0' . $DIGITS . 'u';
        $blk_no       = sprintf($sf_format, $blk_no) or die "Can't zero-pad block number $blk_no\n";

        # Make safely-named block FASTA output file:
        my $seq_file  = $output_name . q{.} . $blk_no . '.fa';
        $seq_file     = safename($seq_file);
        open my $SEQ_ID, '>', $seq_file or die "Can't open sequence output file $seq_file: $!";

        # Print sequences to their appropriate block files:
        foreach my $seq (@seq_names) { 
            if (! exists $data_ref->{'seq_name'}->{$seq}->{'block'} ) { 
                die "Can't identify which block into which to print $seq\n";
            }
            if ( $data_ref->{'seq_name'}->{$seq}->{'block'} == $block ) { 
                print $SEQ_ID '>', $data_ref->{'seq_name'}->{$seq}->{'header'}, "\n";
                my $curr_seq = $data_ref->{'seq_name'}->{$seq}->{'seq'};
                @output_lines = unpack( "a60" x ( length($curr_seq)/60 + 1), $curr_seq );
                foreach $output_line (@output_lines) {
                    print $SEQ_ID "$output_line\n" if ($output_line =~ /\S/);
                }
                # Once all that is done, mark the record as printed (so we can check for unprinted data later):
                $data_ref->{'seq_name'}->{$seq}->{'printed'} = 1;
            }
        }
        close $SEQ_ID or die "Can't close filehandle to sequence output file $seq_file: $!";
    }

    # If *any* data remains unprinted, loudly die:
    my @unprinted_seqs = ();
    foreach my $seq (@seq_names) { 
        if (! exists $data_ref->{'seq_name'}->{$seq}->{'printed'} ) {
            push @unprinted_seqs, $seq;
        }
    }
    if (@unprinted_seqs) { 
         die "Failed to print the following sequences: @unprinted_seqs\n";
    }

    # If all data has been properly printed, then quietly end the successful job:
    exit;
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

sub check_count_val {
    my $val  = $_[0];
    my $type = $_[1];
    if (! looks_like_number($val) ) {
        die "The $type count $val does not look like a number\n";
    }
    if ( ( $val != int($val) ) or ( $val < 2 ) ) {
        die "The $type count $val should give a positive integer, which should be at least 2\n";
    }
    return;
}

