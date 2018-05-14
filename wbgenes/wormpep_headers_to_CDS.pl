#!/usr/bin/perl -w

# Program: wormpep_headers_to_CDS.pl
#    Erich Schwarz, 2/01/02
# 
# Purpose: Get CDS list out of wormpep53 headers.

# If not given a wormpep file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input wormpepN file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".CDS_list");
print "The input file is $infile; the output file $outfile\n";

# Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found. $!\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname. $!\n";

# Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
	if ($_ =~ /^>(\S+)\s+/)
	{
		print OUTFILE "$1\n";
	}
}

close INFILE;
close OUTFILE;
