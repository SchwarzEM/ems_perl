#!/usr/bin/env perl

# censor_small_contigs.pl -- Erich Schwarz <ems394@cornell.edu>, 4/26/2015.
# Purpose: given a WGS, censor tiny contigs (user-specified lower threshold for NON-'tiny'); by default, only remove the 3' terminal contig but leave the 5' terminal contig inviolate, so that annotations of genes won't (generally) be wrecked; however, do allow user to trim off 5' terminal contigs if specifically requested.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles         = ();
my @input_scaffolds = ();
my $scaffold        = q{};
my $threshold       = 0;
my $five_prime;
my $polish;

my $data_ref;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'threshold=i'  => \$threshold,
             'five_prime'   => \$five_prime,
             'polish'       => \$polish,
             'help'         => \$help,   );

if ( $help or (! @infiles) or (! looks_like_number($threshold) ) or ( $threshold != int $threshold ) or ( $threshold <= 0 ) ) { 
    die "Format: censor_small_contigs.pl\n",
        "    --infile|-i       <input stream/files>\n",
        "    --threshold|-t    [minimum length of ACGTacgt contigs to keep, i.e., NOT mask as N; must be positive integer]\n",
        "    --five_prime|-f   [as non-default behavior, trim 5'-terminal small contigs as well as 3'-terminal small contigs]\n",
        "    --polish|-p       [as non-default behavior, trim 5'- and 3'-terminal N/n residues after masking is done]\n",
        "    --help|-h         [print this message]\n",
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

    LOOP: foreach my $scaf (@input_scaffolds) {
        my $orig_seq        = $data_ref->{'scaffold'}->{$scaf}->{'seq'};
        my $orig_header     = $data_ref->{'scaffold'}->{$scaf}->{'header'};
        my $new_seq         = q{};
        my $seen_first_acgt = q{};

        # If there is somehow trailing 3' end, trim it and emit a warning message.
        # Trim trailing 3' end Nn residues before deciding if something's trivial.

        if ( $orig_seq =~ /([Nn]+) \z/xms ) {
            my $trailing_n_residues = $1;
            my $length_n_res = length($trailing_n_residues);
            $orig_seq =~ s/[Nn]+\z//;
            warn "Deleted $length_n_res N/n residues from the 3'-terminus of the scaffold $scaf\n";
        }

        # Censor loudly if $orig_seq somehow *starts* with N; that's just pathological.

        if ( $orig_seq =~ /\A [Nn] /xms ) { 
            warn "Not processing pathological scaffold $scaf, because it starts with one or more N/n residues\n";
            next LOOP;
        }

        # Trivial situation.  Return $new_seq identical to $orig_seq.

        elsif ( $orig_seq =~ /\A [^Nn]+ \z/xms ) { 
            $new_seq = $orig_seq;
        }
 
        # Non-trivial situation.
        # Split into array of very first ACGTacgt, N, etc., last ACGTacgt.
        # End up with $new_seq in which small contigs are masked with N.

        else {
            my @seq_blocks = ();
            if ( $orig_seq !~ /\A [^Nn]+ [Nn]+ [^Nn]+ .* \z/xms ) {
                die "Failed to parse sequence of $scaf\n";
            }
            while ( $orig_seq =~ /\A ([^Nn]+) ([Nn]+) ([^Nn]+ .*) \z/xms ) {
                my $first_acgt = $1;
                my $n_scaffold = $2;
                my $residuum   = $3;

                my $acgt_len = length($first_acgt);

                # Unless a protected 5'-most contig, turn ACGTacgt into N for really small contigs.
                if ( ( $acgt_len < $threshold ) and ( $five_prime or $seen_first_acgt ) ) {                             
                    $first_acgt = ('N' x $acgt_len);
                }

                # Always do this, so only the 5'-most contig is possibly protected.
                $seen_first_acgt = $first_acgt;

                # Process existing sequence blocks:
                push @seq_blocks, $first_acgt, $n_scaffold;
                $orig_seq = $residuum;
            }

            # After that is all done, check the last bit of sequence for being a too-small contig, N-mask if necessary, and then add it in.

            my $residual_len = length($orig_seq);
            if ( $residual_len < $threshold ) {
                $orig_seq = ('N' x $residual_len);
            }
            push @seq_blocks, $orig_seq;

            # Created new single sequence by rejoining the array of sequence blocks.
            $new_seq = join q{}, @seq_blocks;

            # Retrim the 5' and 3' ends if they are now trailing N, and we have specifically asked to polish them off.
            if ($polish) {
                $new_seq =~ s/\A[Nn]+//;
                $new_seq =~ s/[Nn]+\z//;
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
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

