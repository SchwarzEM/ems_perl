#!/usr/bin/env perl

# fasta_size_subset.pl -- Erich Schwarz <ems@emstech.org>, 2/10/2014.
# Purpose: either censor seqs. that are too short (e.g., zero-length) or too long, or, give fraction adding up to a particular size (from largest or smallest seqs.).

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @input_files = ();

my $sequence_size;
my $fraction_size;

my $min;
my $max;

my $output_line     = q{};
my @output_lines    = (); 
my $seq_name        = q{};
my %seqs2headers    = (); 
my %sequences       = ();

my @chosen_seqs = ();

my $prefix = q{};
my $suffix = q{};

my %suf2val = ( K => 1e+03,
                M => 1e+06,
                G => 1e+09, );

my $help;

# To be used as filehandle later:
my $INPUT_FILE;

GetOptions ( 'input_files=s{,}' => \@input_files,
             'sequence_size'    => \$sequence_size,
             'fraction_size'    => \$fraction_size,
             'min=s'            => \$min,
             'max=s'            => \$max, 
             'help'             => \$help, );

if (    $help 
     or (! @input_files                           )
     or ( (! $sequence_size) and (! $fraction_size) ) 
     or ( $sequence_size and $fraction_size         ) ) { 
    die "\n",
        "Format: fasta_size_subset.pl\n",
        "        --input_files|-i    [input file(s), or '-' if stream]\n",
        "        --sequence_size|-s  [boolean: give all sequences (regardless of what they add up to) less than or equal to --max nt, or greater than or equal to --min nt, in size]\n",
        "          [or]\n",
        "        --fraction_size|-f  [boolean: for a given size of of X nt (e.g., 315M), give sequences adding up to that size, selected from either the largest or the smallest sequences:\n",
        "                                --max 'X' gives X nt, built from the largest sequences descending in size;\n",
        "                                --min 'X' gives X nt, built from the smallest sequences ascending in size;\n",
        "\n",
        "        --min [for --sequence_size, non-negative integer]\n",
        "        --max [for --sequence_size, positive integer]\n",
        "              [for --fraction_size, --min or --max must be positive integer/float, optionally with K, M, or G suffix]\n",
        "        --help|-h\n",
        "\n",
        ;
}

# Accept either a stream from '-' or a standard file.
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
            $seq_name = $2;
            $seqs2headers{$seq_name} = $1;
            $sequences{$seq_name} = q{};
        }
        elsif ( $input_line =~ /\S/xms ) {
            $sequences{$seq_name} .= $input_line;
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!\n";
}

# After having imported all sequence data, parse it out 
#     either by min/max contig size, 
#     or by size fractions made from the smallest possible ('min') 
#         or the largest possible ('max') contigs in the input set.

if ($sequence_size) {
    if ( (defined $min) and ( ( $min < 0 ) or ( $min != int($min) ) ) ) { 
        die "Minimum value $min should be a nonnegative integer!\n";
    }
    if ( (defined $max) and ( ( $max < 1 ) or ( $max != int($max) ) ) ) {
        die "Maximum value $max should be a positive integer!\n";
    }
    if ( (defined $min) and (defined $max) and ($min > $max) ) { 
        die "Minimum value $min should be less than or equal to",
            " maximum value $max!\n",
            ;
    }
    # Weed out bad sequences; print out the good ones.
    foreach my $seq_name2 (sort keys %sequences) { 
        my $length = length ( $sequences{$seq_name2} );
        if (  (  (! defined $min ) 
                 or 
                 ( ( defined $min ) and ( $length >= $min ) )  ) 
              and 
              (  (! defined $max )
                 or
                 ( ( defined $max ) and ( $length <= $max ) )  ) 
               ) {
            push @chosen_seqs, $seq_name2;
        }
    }
}

if ($fraction_size) { 
    my $frac_size = 0;
    my $sum_sizes = 0;

    # Do sanity checks, and convert K/M/G abbreviations into real numbers.
    if ( ( defined $min ) and ( defined $max ) ) { 
        die "For the size-fraction option, only one --min or --max fraction argument works.\n";
    }
    if ( defined $min ) {  
        $frac_size = $min;
    }
    if ( defined $max ) { 
        $frac_size = $max;
    }
    if (! looks_like_number($frac_size) ) {
        if ( $frac_size =~ /\A (\S+) (K|M|G) \z /xms ) { 
            $prefix = $1;
            $suffix = $2;
            $frac_size = ( $prefix * $suf2val{$suffix} );
        }
        if (! looks_like_number($frac_size) ) {
            die "Can't parse fraction size argument $frac_size!\n";
        }
    }
    if ( ( $frac_size <= 0 ) or ( $frac_size != int($frac_size) ) ) {
            die "Size-fraction is $frac_size, but instead must be a positive integer of nt.\n";
    }

    # At last, do something with all this...

    my @seq_names = sort { length($sequences{$b}) <=> length($sequences{$a}) } keys %sequences;
    if (defined $min) {
        @seq_names = reverse @seq_names;
    }
    ADD_SEQS: 
    foreach my $seq_name3 (@seq_names) { 
        if ( $sum_sizes < $frac_size ) { 
            push @chosen_seqs, $seq_name3;
            $sum_sizes += length($sequences{$seq_name3});
        }
        else { 
            last ADD_SEQS;
        }
    }
}

# Given a list from somewhere, if it has members, print it out:

foreach my $chosen_seq (@chosen_seqs) { 
    print ">$seqs2headers{$chosen_seq}\n";
    @output_lines = unpack( "a60" x (length($sequences{$chosen_seq})/60 + 1), $sequences{$chosen_seq} );
    foreach $output_line (@output_lines) {
        if ($output_line =~ /\S/) {
            print "$output_line\n";
        }
    }
}

