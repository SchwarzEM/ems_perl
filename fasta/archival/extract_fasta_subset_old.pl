#!/usr/bin/perl -w

# extract_fasta_subset.pl
# Erich Schwarz, 6/26/02

# Purpose: given a list of sequences, extract that set from a bigger FASTA file.

print "List of sequences to extract as FASTA subset? ";
$input_list = <STDIN>;
chomp($input_list);

print "FASTA file to extract subset from? ";
$input_fasta = <STDIN>;
chomp($input_fasta);

$output_subset = $input_fasta . ".subset";

open (INPUT_LIST, "$input_list") || die "Can't open $input_list. $!\n";
open (OUTPUT_SUBSET, ">$output_subset") || die "Can't open $output_subset. $!\n";

while (<INPUT_LIST>) 
{
     $input_name = $_;
     chomp($input_name);
     open (INPUT_FASTA, "$input_fasta") || die "Can't open $input_fasta. $!\n";
     $reading_subset = "no";
     $warn_seq_missing = "yes";
     $warn_seqbody_missing = "yes";
     while (<INPUT_FASTA>) 
     {
         $fasta_line = $_;
         if ($fasta_line =~ /^>$input_name\s/) 
         {
             print OUTPUT_SUBSET "$fasta_line";
             $reading_subset = "yes";
             $warn_seq_missing = "no";
         }
         elsif ($fasta_line =~ /^>\S+/) 
         {
             $reading_subset = "no";
         }
         elsif ($reading_subset eq "yes") 
         {
             print OUTPUT_SUBSET "$fasta_line";
             $warn_seqbody_missing = "no";
         }
     }
     close INPUT_FASTA;
     if ($warn_seq_missing eq "yes") 
     {
         print "Header for sequence $input_name not found in $input_fasta file.\n";
     }
     if ($warn_seqbody_missing eq "yes" && $warn_seq_missing eq "no") 
     {
         print "Header present, but sequence body absent, for\n";
         print "sequence $input_name in $input_fasta file.\n";
     }
}

close INPUT_LIST;
close OUTPUT_SUBSET;
