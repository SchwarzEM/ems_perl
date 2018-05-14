#!/usr/bin/perl -w

# extract_flatfile.pl
# Erich Schwarz, 2/9/04

# Purpose: given a list of flat file entries, extract that set from a flatfile.

chomp($date = `date +"%s"`);

print "List of entries to extract as flatfile subset? ";
chomp($input_list = <STDIN>);

print "Flatfile to extract subset from? ";
chomp($input_flatfile = <STDIN>);

print "Standard text (e.g., '>', 'Locus : ')\n";
print "  that shows the start of a new flat file entry? ";
chomp($start_marker = <STDIN>);

$output_subset = $input_flatfile . $input_list . "." . $date . ".subset";

open (INPUT_LIST, "$input_list") || die "Can't open $input_list. $!\n";
open (OUTPUT_SUBSET, ">$output_subset") || die "Can't open $output_subset. $!\n";

while (<INPUT_LIST>) 
{
     $input_name = $_;
     chomp($input_name);
     open (INPUT_FLATFILE, "$input_flatfile") || die "Can't open $input_flatfile. $!\n";
     $reading_subset = "no";
     $warn_seq_missing = "yes";
     $warn_seqbody_missing = "yes";
     while (<INPUT_FLATFILE>) 
     {
         $flatfile_line = $_;
         if ($flatfile_line =~ /^\W*$input_name\W*/) 
         {
             print OUTPUT_SUBSET "$flatfile_line";
             $reading_subset = "yes";
             $warn_seq_missing = "no";
         }
         elsif ($flatfile_line =~ /^$start_marker/) 
         {
             $reading_subset = "no";
         }
         elsif ($reading_subset eq "yes") 
         {
             print OUTPUT_SUBSET "$flatfile_line";
             $warn_seqbody_missing = "no";
         }
     }
     close INPUT_FLATFILE;
     if ($warn_seq_missing eq "yes") 
     {
         print "Header for sequence $input_name not found in $input_flatfile file.\n";
     }
     if ($warn_seqbody_missing eq "yes" && $warn_seq_missing eq "no") 
     {
         print "Header present, but sequence body absent, for\n";
         print "sequence $input_name in $input_flatfile file.\n";
     }
}

close INPUT_LIST;
close OUTPUT_SUBSET;
