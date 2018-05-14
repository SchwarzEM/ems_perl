#!/usr/bin/perl -w

# Program: wormpep53_ABC_headers_to_loci.pl
#    Erich Schwarz, 6/29/01
# 
# Purpose: Get locus list out of wormpep53 headers.
# Note: this version selects *only* CDSes that end in "A,B,C..."
# to help check against redundant locus names (versus non ABC).

# 1. If not given a wormpep file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input GenPept file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".ABC_locus_list");
print "The input file is $infile; the output file $outfile\n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
	if ($_ =~ /^>\w+\.\d+[A-Z]+ .+ locus:(\S+)\s+/) 
	{
		print OUTFILE "$1\n";
	}
}
