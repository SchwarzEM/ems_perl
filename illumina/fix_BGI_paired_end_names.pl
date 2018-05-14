#!/usr/bin/env perl

# fix_BGI_paired_end_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/31/2011.
# Purpose: give BGI-style FASTQ files '#0/1', '#0/2' as suffixes of paired-end reads, so that bowtie can actually handle them, e.g., for filtering against an unwanted contaminant or background sequence.

use strict;
use warnings;
use Getopt::Long;

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts   = ();

GetOptions ( 'fastq=s' => \$opts{'input_fastq'}, 
             'help'    => \$opts{'help'},         );

if ( $opts{'help'} or (! $opts{'input_fastq'}) ) { 
    die "\n",
        "Format: fix_BGI_paired_end_names.pl\n",
        "          --fastq|-f    [large BGI-style FASTQ to rename so that bowtie doesn't choke on the paired-end read names (opt. '-')]\n",
        "          --help\n",
        "\n",
        ;
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FASTQ;
if ($opts{'input_fastq'} eq '-') {
    # Special case: get the stdin handle
    $INPUT_FASTQ = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FASTQ, '<', $opts{'input_fastq'} or die "Can't open $opts{'input_fastq'}. $!\n";
}

while (<$INPUT_FASTQ>) {
    $head1 = $_;
    $nt_seq = <$INPUT_FASTQ>;
    $head2  = <$INPUT_FASTQ>;
    $quals  = <$INPUT_FASTQ>;

    my $front_header  = q{};
    my $mid_header    = q{};
    my $serial_number = q{};

    chomp( $head1, $nt_seq, $head2, $quals );
    if ( ( $head1 !~ /\A @ \S+ \# [ACGTN]+ \/ [12] \s* \z /xms ) or ( $head2 !~ /\A \+ \s* \z /xms ) ) { 
        warn "Can't parse one or both of these headers:\n";
        warn "$head1\n";
        warn "$head2\n";
        die;
    }
    if ( $head1 =~ /\A (@ \S+) \# ([ACGTN]+) (\/ [12]) \s* \z /xms ) { 
        $front_header  = $1;
        $mid_header    = $2;
        $serial_number = $3;
        $head1 = $front_header . q{_} . $mid_header . q{#0} . $serial_number;
    }
    print "$head1\n", "$nt_seq\n", "$head2\n", "$quals\n", ;
}
close $INPUT_FASTQ or die "Can't close filehandle to FASTA file $opts{'input_fastq'}: $!";

