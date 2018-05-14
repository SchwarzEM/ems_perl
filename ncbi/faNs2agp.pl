#!/usr/bin/env perl

# faNs2agp.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/21/2010.
# Purpose: convert >= 1 FASTA files (with >=1 seqs.) with N runs to AGP, using scaffold_\d+ contig names as in "faNs2ncbi.pl -c" output.

use strict;
use warnings;
use Getopt::Long;

my $scaffold_name = q{};
my $contig_name   = q{};
my %sequences     = ();

my @infiles;
my $minimum_Nrun;
my $help;

my ( $object_name,        $object_start_nt,  $object_end_nt,
     $object_annot_line,  $component_type,   $component_name,
     $component_start_nt, $component_end_nt, $orientation,    );

GetOptions ( 'help'          => \$help,
             "infiles=s{,}"  => \@infiles,
             'minimum=i'     => \$minimum_Nrun, );

if ( $help or (! @infiles ) ) { 
    warn "Format: faNs2agp.pl\n";
    warn "   -i|--infiles  [1+ input files]\n";
    warn "   -m|--minimum  [integer, minimum Ns that split to contigs/split-FASTA: default 10; at least 1]\n";
    die  "   -h|--help\n";
}

# Enforce default value to minimum run of Ns required to split a sequence into contigs.
# The minimum value must be a defined positive integer; otherwise it defaults to 10.
$minimum_Nrun ||= 10;
if ( ( $minimum_Nrun != int($minimum_Nrun) ) or ( $minimum_Nrun < 1 ) ) { 
    warn "Setting --minimum Nrun value to default of 10\n";
    $minimum_Nrun = 10;
}

foreach my $infile (sort @infiles) { 
    open my $INFILE, '<', $infile or die "Can't open input: $infile\n";
    while (my $input_line = <$INFILE>) { 
        chomp $input_line;
        if ($input_line =~ /\A > (\S+) .* /xms) { 
            $scaffold_name = $1; 
            $sequences{$scaffold_name} = q{};
        }
        elsif ( $input_line =~ /\A > /xms ) { 
            die "Aberrant input line: \"$input_line\"\n";
        }
        elsif ( $input_line =~ /[a-zA-Z]/xms ) { 
            $sequences{$scaffold_name} .= $input_line;
        }
    }
    close $INFILE or die "Can't close filehandle to input: $infile\n";
}

foreach my $seq_name2 (sort keys %sequences) { 
    # Do all editing and error-checking with final raw sequence -- this makes catching internal stop codons much easier!
    # Trim terminal stop codons.
    $sequences{$seq_name2} =~ s/\*\z//;

    # Die loudly if given an internal stop codon.
    if ( $sequences{$seq_name2} =~ /\*/xms ) { 
        die "$seq_name2 contains an internal stop codon!\n";
    }

    # Remove all whitespace.
    $sequences{$seq_name2} =~ s/\s//g;

    # Die loudly if given any other weird characters in the sequence.
    my $badchar = '[unclear]';
    if ( $sequences{$seq_name2} =~ /([^a-zA-Z])/xms ) { 
        $badchar = $1;
        die "$seq_name2 contains at least one unacceptable character: \"$badchar\"!\n"; 
    }

    my @contig_lengths   = map { length($_) } split /[nN]{$minimum_Nrun,}/, $sequences{$seq_name2};
    my @contig_names     = ();
    my $contig_number    = @contig_lengths;
    my $contig_digits    = length($contig_number);
    my $contig_length    = 0;

    # Tracked values, for working through multicomponent objects:
    my $current_start_nt = 0;
    my $current_end_nt   = 0;
    my $annot_line       = 0;

    if ( $contig_number == 1 ) {
        # If only one contig in scaffold, life is simple; just annotate the scaffold as a contig.
        $object_name        = $seq_name2;
        $object_start_nt    = 1;
        $object_end_nt      = length($sequences{$seq_name2});
        $object_annot_line  = 1;
        $component_type     ='W';
        $component_name     = $seq_name2;
        $component_start_nt = 1;
        $component_end_nt   = $object_end_nt;
        $orientation        = '+';

        # Format and print contig annotation line:
        print_annot_line();
    }

    # If >= 2 contigs, make a series of annotation lines for contigs and gaps.
    if ( $contig_number >= 2 ) {

        # Make a series of contig names, distinct from the scaffold name.
        foreach my $i (1..$contig_number) {
            my $format = '%0' . $contig_digits . 'u';
            my $contig_tag  = sprintf "$format", $i;
            $contig_name = $seq_name2 . '_' . $contig_tag;
            push @contig_names, $contig_name;
        }

        # Demand that the parsing work.
        if ( $sequences{$seq_name2} !~ / \A [^nN]+ (?: [nN]{$minimum_Nrun,} [^nN]+ )+ \z /xms ) {
            die "Cannot parse sequence of $seq_name2!\n";
        }

        # Get the gap lengths (we already have the contigs themselves, in @contigs).
        # Use 'while ... /g' to do progressive matches through the string.

        my $gap_Ns      = q{};
        my @gap_lengths = ();
        while ( $sequences{$seq_name2} =~ / ([nN]{$minimum_Nrun,}) [^nN]+ /gxms ) {
            $gap_Ns = $1;
            push @gap_lengths, length($gap_Ns);
        }

        # Now that we have contig names, and contig and gap lengths, produce annotation lines.
        foreach my $i (1..$contig_number) {

            # Use for Perl indexes, notably in gap contig lengths list.
            my $j = ($i - 1);

            # Start off the object-oriented coordinates and line number.
            if ( $i == 1 ) {
                $current_start_nt  = 1;
                $annot_line        = 1;
            }

            # Now that we know the next length, increment $current_end_nt up.
            $current_end_nt += $contig_lengths[$j];

            # Write up contig annotation:
            $object_name        = $seq_name2;
            $object_start_nt    = $current_start_nt;
            $object_end_nt      = $current_end_nt;
            $object_annot_line  = $annot_line;
            $component_type     ='W';
            $component_name     = $contig_names[$j];    # Not $i!
            $component_start_nt = 1;
            $component_end_nt   = $contig_lengths[$j];
            $orientation        = '+';
         
            # Format and print contig annotation line:
            print_annot_line();

            # Immediately increment $current_start_nt and $annot_line.
            # (Leave $current_end_nt unchanged until we know the next length.)
            $current_start_nt  += $contig_lengths[$j];
            $annot_line++;

            # If there's a following gap, annotate it too:
            if ( $i < $contig_number ) { 
                # Update:
                $current_end_nt += $gap_lengths[$j];

                # Write up gap annotation:
                $object_name        = $seq_name2;
                $object_start_nt    = $current_start_nt;
                $object_end_nt      = $current_end_nt;
                $object_annot_line  = $annot_line;
                $component_type     = 'N';
                $component_name     = $gap_lengths[$j];
                $component_start_nt = 'fragment';
                $component_end_nt   = 'yes';
                $orientation        = q{};

                # Format and print gap annotation line:
                print_annot_line();

                # Update:
                $current_start_nt  += $gap_lengths[$j];
                $annot_line++;
            }
        }
    }
}

sub print_annot_line { 
    my $_annotation  =   $object_name
                       . "\t"
                       . $object_start_nt
                       . "\t"
                       . $object_end_nt
                       . "\t"
                       . $object_annot_line
                       . "\t" 
                       . $component_type
                       . "\t" 
                       . $component_name
                       . "\t" 
                       . $component_start_nt
                       . "\t" 
                       . $component_end_nt
                       . "\t" 
                       . $orientation
                       ;

    print "$_annotation\n";
    return;
}

