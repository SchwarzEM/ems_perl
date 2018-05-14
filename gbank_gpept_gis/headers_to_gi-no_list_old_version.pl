#!/usr/bin/perl -w

# Program: headers_to_gi-no_list.pl
#    Erich Schwarz, 3/12/01
# 
# Purpose: To take the hit list from BlastP and make it a useful
#    gi number list for batch Entrez.
#
#    (Whose GenPept downloads can then be processed with 
#     genpept_to_clean_tfa.pl.)

# 1a. If not given a BlastP report, or excerpt of same, as an argument, 
#         ask for its name.

if ($#ARGV != 0) {
      print "Required: input BlastP results file!\n";
      print "What will input file be? ";
      $infile = <STDIN>;
} else {
      $infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".gi-nos");
print "The input file is $infile; the output file $outfile\n";

# Initialize a control variable which I'll use to prevent
# text from being scanned as coding sequences.

$aa_read_toggle = 0;

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "GenPept file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
    if ($_ =~ /^gi\|([0-9]*)\|/) 
    {
        print OUTFILE ($1);
        print OUTFILE ("\n");
    }
}
