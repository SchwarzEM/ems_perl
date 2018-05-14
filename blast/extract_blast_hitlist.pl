#!/usr/bin/perl -w

# extract_blast_hitlist.pl
# Erich Schwarz, 5/14/02.
# Designed to extract useful hitlists out of Blast searches
#   (e.g., mass BlastX searches with ESTs versus various databases).

# Datestamping is used to ensure that a given analyses' files
#   are kept sorted out; two different runs of the script will
#   automatically generate distinctly stamped files and directories.

print "\n";
print "This is most useful on mass Blast searches, e.g., mass\n";
print "BlastX searches with ESTs versus various databases.\n";
print "\n";
print "Blast search to extract?:   ";

$input_file=<STDIN>;
chomp($input_file);

$list_of_hits = $input_file . ".list_of_hits";

open (HITLIST, ">$list_of_hits");

# [Format of things to extract:]
# Database: ../nr
#            930,802 sequences; 292,018,088 total letters
# 
# Searching..................................................done
# 
#                                                                    Score     E
# Sequences producing significant alignments:                        (bits)  Value
# 
# gi|18543885|ref|XP_086920.1| (XM_086920) similar to putative [Ho...    78  1e-14
# gi|18543890|ref|XP_086919.1| (XM_086919) similar to ribosomal pr...    78  1e-14
# gi|18543881|ref|XP_086921.1| (XM_086921) similar to ribosomal pr...    78  2e-14
# gi|18595696|ref|XP_060067.2| (XM_060067) similar to ribosomal pr...    67  3e-10
# 
# >gi|18543885|ref|XP_086920.1| (XM_086920) similar to putative [Homo sapiens]
#           Length = 232


$scan_hits = "no";
$print_next_line = "no";

open (INPUT_FILE, "$input_file");

while (<INPUT_FILE>) 
{
    $input_line = $_;
    chomp($input_line);

    if ($input_line =~ /^Sequences producing significant alignments/)
    {
        $scan_hits = "yes";
    }
    elsif ($input_line =~ /Query= /)
    {
        $print_next_line = "yes";
        print HITLIST "\n";
        print HITLIST "$input_line";
    }
    elsif ($print_next_line eq "yes")
    {
        $print_next_line = "no";
        print HITLIST " ";
        print HITLIST "$input_line\n";
    }
    elsif ($input_line =~ /^Database: /) 
    {
        print HITLIST "$input_line\n";
    }
    elsif ($input_line =~ /No hits found/)
    {
        print HITLIST "[ no hits found ]\n";
    }
    elsif ($input_line =~ /^>\S+/) 
    {
        $scan_hits = "no";
    }
    elsif ($input_line =~ /^gi\|\d+\|/) 
    {
        if ($scan_hits eq "yes") 
        {
            print HITLIST "$input_line\n";
        }
    }
    elsif ($input_line =~ /^\S+\s+/) 
    {
        if ($scan_hits eq "yes") 
        {
            print HITLIST "$input_line\n";
        }
    }
}

close HITLIST;

print "\n";
print "Job done.\n";
print "\n";
