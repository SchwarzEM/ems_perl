#!/usr/bin/perl -w

# Program: wormpepNN_nonABC_headers_to_loci.pl
#    Erich Schwarz, 8/29/01
# 
# Purpose: Get locus list out of wormpepNN headers.
# Note: this version selects against CDSes that end in "A,B,C...",
# in order to check against redundant locus names.

# 1. If not given a wormpep file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input wormpepNN file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".nonABC_locus_list");
print "The input file is $infile; the output file $outfile\n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
	if ($_ =~ /^>\w+\.\d+ .+ locus:(\S+)\s+/) 
	{
		print OUTFILE "$1\n";
	}
}
