#!/usr/bin/perl -w

# Program: headers_to_gi-no_list.pl
#    Erich Schwarz, 11/25/01
# 
# Purpose: To take either a hitlist from BlastP or gi-numbered FASTA 
#    file headers, and extract from either one a useful
#    gi number list for batch Entrez.
#
#    (Whose GenPept downloads can then be processed with 
#     genpept_to_clean_tfa.pl.)

if ($#ARGV != 0) 
{
    print "What will input file be? ";
    $infile = <STDIN>;
} 
else 
{
    $infile = $ARGV[0];
}

chomp ($infile);
$rough_outfile = ($infile . "rough_outfile");
$outfile = ($infile . ".gi-nos");

open (INFILE, "$infile") || die "File $infile not found. $!\n";
open (ROUGH_OUTFILE, ">$rough_outfile") || die "Couldn't open working file $rough_outfile. $!\n";

while (<INFILE>) 
{
    if ($_ =~ /^gi\|([0-9]+)\|/) 
    {
        print ROUGH_OUTFILE ($1);
        print ROUGH_OUTFILE ("\n");
    }
    elsif ($_ =~ /^>([0-9]+)[^0-9]+/) 
    {
        print ROUGH_OUTFILE ($1);
        print ROUGH_OUTFILE ("\n");
    }
}
close INFILE;
close ROUGH_OUTFILE;

system "sort $rough_outfile | uniq > $outfile";
system "rm $rough_outfile";

print "The input file is $infile; the sort-ed and uniq-ed output file is $outfile.\n";
