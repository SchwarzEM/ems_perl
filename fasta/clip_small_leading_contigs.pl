#!/usr/bin/env perl

# clip_small_leading_contigs.pl -- Erich Schwarz <ems394@cornell.edu>, 7/22/2013.
# Purpose: given a genome assembly with some remaining leading very small contigs, clip scaffolds to remove them, and report the nt shifted for each scaffold.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles         = ();
my @input_scaffolds = ();
my $scaffold        = q{};
my $threshold       = 0;
my $extra_nt        = 0;

my $data_ref;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'threshold=i'  => \$threshold,
             'help'         => \$help,   );

if ( $help or (! @infiles) or (! looks_like_number($threshold) ) or ( $threshold != int $threshold ) or ( $threshold <= 0 ) ) { 
    die "Format: clip_small_leading_contigs\n",
        "    --infile|-i     <input stream/files>\n",
        "    --threshold|-t  [maximum length of leading 5' ACGTacgt contigs to censor (along with first block of scaffolding N; must be positive integer]\n",
        "    --help|-h       [print this message]\n",
        ;
}

foreach my $infile (@infiles) { 
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( $input =~ /\A > (\S+) .*\z/xms ) { 
            $scaffold = $1;
            if ( exists $data_ref->{'scaffold'}->{$scaffold} ) { 
                die "Redundant sequence name: $scaffold\n";
            }
            push @input_scaffolds, $scaffold;

            # Note that $input includes the starting '>' for a FASTA record line.
            $data_ref->{'scaffold'}->{$scaffold}->{'header'} = $input;
        }
        elsif ( $input =~ /\A > /xms ) {
            die "Can't parse input line: $input\n";
        }
        else {
            $data_ref->{'scaffold'}->{$scaffold}->{'seq'} .= $input;
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

LOOP: foreach my $scaf (@input_scaffolds) {
    my $orig_seq    = $data_ref->{'scaffold'}->{$scaf}->{'seq'};
    my $orig_header = $data_ref->{'scaffold'}->{$scaf}->{'header'};
    my $new_seq     = q{};
    my $start_acgt  = q{};

    # Trivial situation.  Return $new_seq identical to $orig_seq.
    if ( $orig_seq =~ /\A [^Nn]+ \z/xms ) { 
        $new_seq = $orig_seq;
    }

    # Non-trivial situation.
    # Split into very first ACGTacgt, first N, etc.
    # End up with $new_seq in which we censor very first ACGTacgt + first N if very first ACGTacgt is too small.

    else {
        # Before we start editing, deal with pathological scaffolds that start with N.
        # For this script, we *do* want to deal with $orig_seq if it starts with N; though it is still a good idea to have a warning message about it. 
        # Since $threshold *must* be at least 'A', prepending 1 nt is guaranteed to make a trimmable 5' end for any scaffold that had previously started with N.
        if ( $orig_seq =~ /\A ([Nn]) /xms ) {
            warn "Pathological scaffold $scaf starts with an N residue; however, it will still be 5'-trimmed of N residues.\n";
            $orig_seq = q{A} . $orig_seq;
            $extra_nt = 1;
        }

        # Next, do the trimming.
        if ( $orig_seq =~ /\A ([^Nn]+) ([Nn]+) ([^Nn]+ .*) \z/xms ) {
            my $first_acgt = $1;
            my $n_scaffold = $2;
            my $residuum   = $3;

            # Default:
            $new_seq = $orig_seq;

            # But if we have a leading small contig, clip, and report the clip:
            my $acgt_len = length($first_acgt);
            my $trim_len = length($n_scaffold);
            $trim_len   += $acgt_len;
            # Unless the scaffold is pathological, $extra_nt is 0; if it is pathological, $extra_nt is 1.
            $trim_len   -= $extra_nt;
            if ( $acgt_len < $threshold ) { 
                $new_seq = $residuum;
                warn "Trimmed 5' end of $scaf by $trim_len nt; correct gene annotations accordingly.\n";
            }
        }
        else { 
            die "Failed to parse sequence of $scaf\n";
        }
    }

    # At last, print the header and sequence.
    # $orig_header kept the '>', remember.
    print "$orig_header\n";
    my @output_lines = unpack("a60" x (length($new_seq)/60 + 1), $new_seq);
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

