#!/usr/bin/perl -w

# Program: wormpep_headers_to_locus.pl
#    Erich Schwarz, 6/29/01
# 
# Purpose: Get locus list out of wormpep53 headers.
# Note: this version selects everybody, including CDSes that 
# end in "A,B,C..." which is not helpful for checking against 
# redundant locus names.

# 1. If not given a wormpep file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input wormpepN file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".locus_list");
print "The input file is $infile; the output file $outfile\n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
	if ($_ =~ /^>.+locus:(\S+)\s+/) 
	{
		print OUTFILE "$1\n";
	}
}
