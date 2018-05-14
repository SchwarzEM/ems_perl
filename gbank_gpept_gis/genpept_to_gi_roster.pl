#!/usr/bin/perl -w

# Program: genpept_to_gi_roster.pl
#    Erich Schwarz, 9/2/99
# 
# Purpose: Starting with a batch GenPept file from batch 
#    Entrez, produce a clean list with one line for each
#    gi number, and with the following on each line:
#
#    gi no.; gene name; product name; species
#
#    This has two uses:
#       
#    In its own right it makes nice summary data from 
#       a batch-Entrez output.
#    The *.gi_roster file produced here can be fed into
#       another script to make tailored FASTA files with
#       optimized (from my perspective) headers.

# 1. If not given a genpept file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input GenPept file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".gi_roster");
print "The input file is $infile; the output file $outfile\n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "GenPept file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
	if ($_ =~ /^PID[^0-9]*([0-9]+)$/) 
	{
		$aa_read_toggle = 0;
		$gi_number = $_;
		chomp($gi_number);
		$gi_number =~ s/[^0-9]//g;
		$species = "no_species_name";
		$gene_name = "no_gene_name";
		$product_name = "no_product_name";
	}
	elsif ($_ =~ /ORGANISM[\s]+[.]*/) 
	{
		$species = $_;
		chomp($species);
		$species =~ s/ORGANISM[\s]+//g;
	}
	elsif ($_ =~ /gene=/) 
	{
		$gene_name = $_;
		chomp($gene_name);
		$gene_name =~ s/[\W]*gene=//;
		$gene_name =~ s/\"//g;
	}
	elsif ($_ =~ /product=/) 
	{
		$product_name = $_;
                chomp($product_name);  
                $product_name =~ s/[\W]*product=//;
                $product_name =~ s/\"//g;
	}
	elsif ($_ =~ /ORIGIN/) 
	{
                print OUTFILE ($gi_number . "  " . $gene_name . "  " . $product_name . "  ");
                print OUTFILE ("gi|" . $gi_number . "  " . $species . "\n");
                $aa_read_toggle = 1;
        }
}
