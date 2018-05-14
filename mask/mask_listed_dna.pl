#!/usr/bin/env perl

# mask_listed_dna.pl -- Erich Schwarz, 11/5/2007; big update on 10/16/2011; <emsch@caltech.edu>.
# Purpose: either use 'X' or 'N' to hardmask listed ranges of nucleotides in one or more FASTA files, or use all-undercase to softmask them.

use strict;
use warnings;
use Getopt::Long;

my $listfile    = q{};
my @fasta_files = ();
my $mask        = q{};
my $inverse;
my $softmask;
my $help;

GetOptions ( 'list=s'     => \$listfile,   
             'fasta=s{,}' => \@fasta_files,
             'mask=s'     => \$mask,
             'inverse'    => \$inverse,
             'soft'       => \$softmask,
             'help'       => \$help,         );

my $seqname             = q{};
my $header              = q{};
my %seqnames2headers    = ();
my $seqnames2res_ref;
my $seqnames2ranges_ref;

if ( $help or ( (! $softmask) and (! $mask) ) or ( $softmask and $mask) or (! $listfile) or (! @fasta_files) ) { 
    die "Format: mask_listed_dna.pl\n",
        "    --list|-l     [list file; for each line: seqname  nt_start  nt_end\n",
        "    --fasta|-f    [FASTA file(s, one or more)]\n",
        "    --mask|-m     [choice of single residue (e.g., X or N) used for hardmasking]\n",
        "    --inverse|-i  [optionally, leave listed residues *unmasked*, masking entire remainder of sequence(s)]\n",
        "    --soft|-s     [use softmasking rather than hardmasking with X|N; also works for inverse]\n",
        "    --help|-h     [print this message]\n",
        ;
}

if (! $softmask) { 
    if ( $mask !~ /\A [A-Za-z] \z /xms ) { 
        die "Mask nucleotide specified as \"$mask\"; it needs to be a single alphabetical letter.\n";
    }
    if ( $mask !~ /\A [XxNn] \z /xms ) {
        warn "Mask nucleotide specified as \"$mask\"; a more standard choice would be \"X\", \"x\", \"N\", or \"n\".\n";
    }
}

foreach my $fasta_input (@fasta_files) { 
    open my $FASTA, '<', $fasta_input or die "Can't open FASTA file $fasta_input: $!";
    while (my $input = <$FASTA>) {
        chomp $input;
        if ($input =~ /\A > ( (\S+) .*) \z /xms) { 
            $header  = $1;
            $seqname = $2;
            if ( exists $seqnames2headers{$seqname} ) { 
                die "Redundant sequence name: $seqname\n";
            }
            $seqnames2headers{$seqname} = $header;
        }
        elsif ($input !~ /\A >/xms) { 
            $input =~ s/\s//g;
            $input =~ s/[^a-zA-Z]//g;

            # For softmasking to work, default must be ALL-CAPS:
            if ($softmask) { 
                $input = uc($input);
            }

            # Then, store as array-ized single residues:
            if ( $input =~ / [a-zA-Z]+ /xms ) { 
                my @residues = split //, $input;
                if (! $seqname) { 
                    die "Trying to list residues of an unnamed sequence!\n";
                }
                push @{ $seqnames2res_ref->{$seqname} }, @residues;
            }
        }
        else { 
            die "Can't parse FASTA line: $input\n";
        }
    }
    close $FASTA or die "Can't close filehandle to FASTA file $fasta_input: $!";
}

# No sequences? game over!
if (! keys %{ $seqnames2res_ref } ) { 
    die "No sequences recorded!\n";
}

open my $LIST, "$listfile" or die "Can't open list file $listfile: $!";
while (my $input = <$LIST>) { 
    chomp $input;

    # Record a range of residues from a named sequence.
    if ( $input =~ /\A (\S+) \s+ (\d+) \s+ (\d+) \s* \z /xms ) { 
        my ($sqnam, $num1, $num2) = ($1, $2, $3);
        # Warn about irrelevant sequences, but don't die:
        if (! exists $seqnames2res_ref->{$sqnam} ) { 
           warn "Sequence $sqnam is in list file but not FASTA file(s)!\n";
        }
        $seqnames2ranges_ref->{$sqnam}->{$num1} = $num2;
    }

    # Enforce clean inputs.
    else {
        die "Malformatted listfile line: $input\n";
    }
}
close $LIST;

foreach my $sqnam (sort keys %{ $seqnames2res_ref }) { 
    # Get length of sequence:
    my $seq_len = @{ $seqnames2res_ref->{$sqnam} };

    # Shift it from natural to Perl 0-based counting:
    $seq_len--;

    # Prepare to mark specific parts of entire sequence:
    my %flagged_residues = ();

    # Flagging residues is optional, not mandatory.
    if ( exists $seqnames2ranges_ref->{$sqnam} ) { 
        # Go through sequence and mask in ascending residue order.
        foreach my $start (sort { $a <=> $b } keys %{ $seqnames2ranges_ref->{$sqnam} } ) { 
            # Define range values.
            my $end = $seqnames2ranges_ref->{$sqnam}->{$start};

            # Shift by -1 to adjust from molbiol to Perl residue-counting.
            $start--;
            $end--;

            # Add to the set of flagged residues, which we will either mask (usually) or spare unmasked (with the 'inverse' option).
            # Note that I am assuming that Perl will cope gracefully with 1-nt sites for which $start equals $end!
            foreach my $i ($start..$end) { 
                $flagged_residues{$i} = 1;
            }
        }
    }

    # But the decision to do something about the sequences isn't optional.
    if ( ( exists $seqnames2ranges_ref->{$sqnam} ) or ( $inverse and (! exists $seqnames2ranges_ref->{$sqnam} ) ) ) { 
        foreach my $seq_residue (0..$seq_len) { 
            if (     ( (! $inverse ) and ( exists $flagged_residues{$seq_residue}  ) ) 
                  or ( $inverse      and (! exists $flagged_residues{$seq_residue} ) ) ) { 
                # Two different options here, which are coded above to never co-occur:
                if ($mask) { 
                    $seqnames2res_ref->{$sqnam}->[$seq_residue] = "$mask";
                }
                if ($softmask) {
                    $seqnames2res_ref->{$sqnam}->[$seq_residue] = lc($seqnames2res_ref->{$sqnam}->[$seq_residue]);
                }
            }
        }
    }
}

# Export results as a single FASTA text.
foreach my $sqnam (sort keys %{ $seqnames2res_ref }) { 
    my $fullseq = join q{}, @{ $seqnames2res_ref->{$sqnam} };
    print ">$seqnames2headers{$sqnam}\n";
    my @output_lines = unpack("a60" x (length($fullseq)/60 + 1), $fullseq);
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/xms) { 
            print "$output_line\n";
        }
    }
}

