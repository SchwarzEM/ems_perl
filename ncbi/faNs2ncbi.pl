#!/usr/bin/env perl

# faNs2ncbi.pl -- Erich Schwarz <ems394@cornell.edu>, 9/20/2010; updated a bit, 8/4/2016.
# Purpose: convert >= 1 FASTA files (with >=1 seqs.) with N runs to either pure-contig FASTA or NCBI's "split FASTA"; also does uniform_fasta.pl reformatting.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $output_line  = q{};
my @output_lines = ();
my $seq_name     = q{};
my %seqs2headers = ();
my %sequences    = ();
my @seq_names    = ();


my $allcaps;
my $contigs;
my @infiles;
my $keep;
my $minimum_Nrun;
my $nopad;
my $split_fasta;

my $help;

GetOptions ( 'allcaps'      => \$allcaps, 
             'contigs'      => \$contigs, 
             'help'         => \$help,
             "infiles=s{,}" => \@infiles,
             'minimum=i'    => \$minimum_Nrun,
             'keep'         => \$keep, 
             'nopad'        => \$nopad,
             'split'        => \$split_fasta,  ); 

if (    $help 
     or (! @infiles                           )
     or ( (! $contigs ) and (! $split_fasta ) ) 
     or ( $contigs and $split_fasta           ) ) { 
    warn "Format: faNs2ncbi.pl\n";
    warn "   -a|--allcaps\n";
    warn "   -c|--contigs  [i.e., pure contigs output] <or> -s|--split [NCBI's split-FASTA output]\n";
    warn "   -i|--infiles  [1+ input files]\n";
    warn "   -m|--minimum  [integer, minimum Ns that split to contigs/split-FASTA: default 10; at least 1]\n";
    warn "   -n|--nopad    [do not pad numbers of contigs with \"0\"]\n";
    warn "   -k|--keep     [keep original order of input sequences]\n";
    die  "   -h|--help\n";
}

# Enforce default value to minimum run of Ns required to split a sequence into contigs.
# The minimum value must be a defined positive integer; otherwise it defaults to 10.
$minimum_Nrun ||= 10;
if ( ( $minimum_Nrun != int($minimum_Nrun) ) or ( $minimum_Nrun < 1 ) ) { 
    warn "Setting --minimum Nrun value to default of 10\n";
    $minimum_Nrun = 10;
}

if (! $keep) {
    @infiles = sort @infiles;
}

foreach my $infile (@infiles) { 
    open my $INFILE, '<', $infile or die "Can't open input: $infile\n";
    while (my $input_line = <$INFILE>) { 
        chomp $input_line;
        if ($input_line =~ /\A > (\S+) (.*) /xms) { 
            $seq_name = $1; 
            $seqs2headers{$seq_name} = $2;
            $sequences{$seq_name} = q{};
            push @seq_names, $seq_name;
        }
        elsif ( $input_line =~ /\A > /xms ) { 
            die "Aberrant input line: \"$input_line\"\n";
        }
        elsif ( $input_line =~ /[a-zA-Z]/xms ) { 
            $sequences{$seq_name} .= $input_line;
        }
    }
    close $INFILE or die "Can't close filehandle to input: $infile\n";
}

if (! $keep) {
    @seq_names = sort keys %sequences;
}

foreach my $seq_name2 (@seq_names) { 
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

    # Optionally, write the sequence in GREAT RUNES.
    if ($allcaps) {
        $sequences{$seq_name2} =~ tr/[a-z]/[A-Z]/;
    }

    my @contigs = split /[nN]{$minimum_Nrun,}/, $sequences{$seq_name2};
    my @contig_names = ();
    my $contig_number = @contigs;
    my $count_length = length($contig_number);

    # Make a list of contig FASTA header names that are appropriate for the desired output.
    # If only one contig in scaffold, life is simple.
    if ( $contig_number == 1 ) {
        push @contig_names, $seq_name2;
    }

    # If >= 2 contigs, different choices depend on the option being used.
    if ( $contig_number >= 2 ) {
        foreach my $i (1..$contig_number) {
            # Two different printing options.
            # For '-c|--contigs-only', prepare to print a string of pure contigs, given minimum size of splitting N run.
            if ($contigs) {
                my $format = '%0' . $count_length . 'u';

                my $contig_tag = $i;
                if (! $nopad) {
                    $contig_tag = sprintf "$format", $i;
                }

                my $contig_name = $seq_name2 . q{.} . $contig_tag;
                push @contig_names, $contig_name;
                $seqs2headers{$contig_name} = $seqs2headers{$seq_name2};
            }
            # For '-s|--split-fasta', prepare to print NCBI's own 'split FASTA' format, given minimum size of splitting N run.
            if ($split_fasta) {
                push @contig_names, $seq_name2;
                my $gap_Ns     = q{};
                my $gap_length = q{};
                my $gap_name   = q{};
                # Use 'while ... /g' to do progressive matches through the string.
                while ( $sequences{$seq_name2} =~ /([nN]{$minimum_Nrun,})[^nN]/gxms ) { 
                    $gap_Ns = $1;
                    $gap_length = length($gap_Ns);
                    $gap_name = '?' . "$gap_length";
                    push @contig_names, $gap_name;
                }
            }
        }
    }

    foreach my $i (1..$contig_number) {
        my $j = $i - 1; 
        print '>';
        print $contig_names[$j];
        print "\t$seqs2headers{$contig_names[$j]}" if ($seqs2headers{$contig_names[$j]});
        print "\n";
        @output_lines 
            = unpack("a60" x (length($contigs[$j])/60 + 1), $contigs[$j]);
        foreach $output_line (@output_lines) { 
            if ($output_line =~ /\S/) { 
                print "$output_line\n";
            }
        }
    }
}

