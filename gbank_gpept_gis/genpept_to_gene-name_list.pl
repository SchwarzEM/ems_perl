#!/usr/bin/perl -w

# Program: genpept_to_gene-name_list.pl
#    Erich Schwarz, 1/21/02
# 
# Purpose: Extract and order a gene name list from GenPept files.

# 1a. If not given a genpept file as argument, ask for its name.

if ($#ARGV != 0) 
{
    print "Required: input GenPept file!\n";
    print "What will input file be? ";
    $infile = <STDIN>;
} 
else 
{
    $infile = $ARGV[0];
}

chomp ($infile);
$rough_genename_outfile = ($infile . ".gene-name-list.1");
$rough_gi_no_outfile = ($infile . ".gi-number-list.1");
$genename_outfile = ($infile . ".gene-name-list");
$gi_no_outfile = ($infile . ".gi-number-list");

print "The input file is $infile;\n";
print "    the output files are $genename_outfile and $gi_no_outfile.\n";
$gene_name = "no_gene_name";
$product_name = "no_product_name";
$gi_number = "";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "Couldn't open GenPept file $infile. $!\n";
open (ROUGH_GENENAME_OUTFILE, ">$rough_genename_outfile") || die "Couldn't open file $rough_genename_outfile. $!\n";
open (ROUGH_GI_NUMBER_OUTFILE, ">$rough_gi_no_outfile") || die "Couldn't open file $rough_gi_no_outfile. $!\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
    if ($_ =~ /^PID[^0-9]*([0-9]+)$/) 
    { 
        $gi_number = $1;
        chomp($gi_number);
        $gi_number =~ s/[^0-9]//g;
        $gene_name = "no_gene_name";
        $product_name = "no_product_name";
    }
    elsif ($_ =~ /gene=/) 
    {
        $gene_name = $_;
        chomp($gene_name);
        $gene_name =~ s/[\W]*gene=//;
        $gene_name =~ s/\"//g;
        $gene_name =~ tr/a-z/A-Z/;
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
         if ($gene_name =~ /\w+\.\S+/) 
         { 
             print ROUGH_GENENAME_OUTFILE "$gene_name\n";
             print ROUGH_GI_NUMBER_OUTFILE "$gi_number\n";
         } 
         elsif ($product_name eq "no_product_name") 
         {
             print ROUGH_GENENAME_OUTFILE "$gi_number\n";
             print ROUGH_GI_NUMBER_OUTFILE "$gi_number\n";
         }
         else 
         {
             print ROUGH_GENENAME_OUTFILE "$product_name\n";
             print ROUGH_GI_NUMBER_OUTFILE "$gi_number\n";
         }
    }
}

close ROUGH_GENENAME_OUTFILE;
close ROUGH_GI_NUMBER_OUTFILE;

system "sort $rough_genename_outfile | uniq - > $genename_outfile";
system "sort $rough_gi_no_outfile | uniq - > $gi_no_outfile";
system "rm $rough_genename_outfile";
system "rm $rough_gi_no_outfile";
