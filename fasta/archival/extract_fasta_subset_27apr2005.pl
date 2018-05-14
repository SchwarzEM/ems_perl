#!/usr/bin/perl

# extract_fasta_subset.pl
# Erich Schwarz, 4/27/05

# Purpose: given a listed set of sequences, extract that set from a bigger FASTA file; warn about missed sequences.

unless ($#ARGV == 1) { 
    die "Format: extract_fasta_subset.pl  [seqs. list]  [large FASTA to extract]\n";
}

chomp ($input_list  = $ARGV[0]);
chomp ($input_fasta = $ARGV[1]);
$output_subset = $input_fasta . "." . $input_list . ".subset";
$warnings      = $input_fasta . "." . $input_list . ".warnings";

open (INPUT_LIST, "$input_list") || die "Can't open $input_list. $!\n";
while (<INPUT_LIST>) {
    chomp;
    $input_names{$_} = 1;
}
close INPUT_LIST;

open (INPUT_FASTA, "$input_fasta") || die "Can't open $input_fasta. $!\n";
while (<INPUT_FASTA>)
{
    $fasta_line = $_;
    if ( ($fasta_line =~ /^>(\S+)\s*/) && ($input_names{$1} == 1) ) {
        print "$fasta_line";
        $reading_subset = "yes";
        $input_names{$1} = 0;
    }
    elsif ( ($fasta_line =~ /^>(\S+)\s*/) && ($input_names{$1} != 1) ) {
        $reading_subset = "no";
    }
    elsif ($reading_subset eq "yes") {
        print "$fasta_line";
    }
}
close INPUT_FASTA ;

# scan hash for any non-zero values, warn that these were missed.

open (WARNINGS, ">$warnings") || die "Can't open $warnings. $!\n";
foreach $key (sort keys %input_names) {
    if ($input_names{$key} == 1) { 
        print WARNINGS "$key not found in $input_fasta\n";
    }
}
