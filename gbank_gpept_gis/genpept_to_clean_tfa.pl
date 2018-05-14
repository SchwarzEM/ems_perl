#!/usr/bin/perl -w

# Program: genpept_to_clean_tfa.pl
#    Erich Schwarz, 9/02/2003
# 
# Purpose: Large-scale homolog searches can be aided by batch Entrez,
#    but neither FASTA nor GenPept downloads really are ideal.
#    The FASTA files give gi|N but not real gene names; they also
#      can have vast comment lines.  The GenPept files have all
#      data I want but not in FASTA format.  Thus, this script is
#      designed to take a GenPept file and convert it into a FASTA
#      file with the following characteristics: the ">" line has
#      the *gene name* as its leading string; the gi|N identifier
#      is given later; and a clean FASTA sequence follows.

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
$outfile = ($infile . ".tfa");
print "The input file is $infile; the output file $outfile\n";

if ($#ARGV != 1) 
{
    print "This can output .tfas named either by their\n";
    print "   gi no. (the default) or accession name.\n";
    print "   Which format is wanted? [\"gi\" or \"acc\"]: ";
    $format = <STDIN>;
} 
else 
{
    $format = $ARGV[1];
}

chomp ($format);

if ($format =~ /^gi/) 
{
    $format_toggle = 0;
    print "The output format is: .tfas named by gi no.\n";
}
elsif ($format =~ /^acc/) 
{
    $format_toggle = 1;
    print "The output format is: .tfas named by accession name.\n";
}
else 
{
    $format_toggle = 0;
    print "The output format BY DEFAULT is: .tfas named by gi no.\n";
}

# Initialize a control variable which I'll use to prevent
# text from being scanned as coding sequences.

$aa_read_toggle = 0;

open (INFILE, "$infile") || die "GenPept file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outfile\n";

while (<INFILE>) 
{
    if ($_ =~ /LOCUS/) 
    {
        $aa_read_toggle = 0;
    }

    elsif ($_ =~ /^ACCESSION\s+(\S+)\s*$/)   # This is the invariant identifier used by NCBI.
    {
        $acc_number = $1;
        $aa_read_toggle = 0;
        chomp($acc_number);
        $species = "no_species_name";
        $gene_name = "no_gene_name";
        $product_name = "no_product_name";
        print OUTFILE (">");
    }

    elsif ($_ =~ /^VERSION.+GI.([0-9]+)$/)   # new format for GI numbers: VERSION     AAF57416.1  GI:7302326
    {
        $gi_number = $1;
        $aa_read_toggle = 0;
        chomp($gi_number);
        $gi_number =~ s/[^0-9]//g;
    }
    elsif ($_ =~ /ORGANISM[\s]+[.]*/) 
    {
        $species = $_;
        chomp($species);
        $species =~ s/ORGANISM[\s]+//g;
        $species =~ s/\<\/a>//;         # this is a very specific kludge, designed in
        $species =~ s/\<a href.+\>//;   # 5/15/02 to remove NCBI cruft from headers...
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
    elsif ($_ =~ /ORIGIN/)  # Direct next output by whether $format_toggle = 1 and whether gene has name.
    {
        if ($format_toggle == 1) 
        {
            print OUTFILE ($acc_number . "   gi|" . $gi_number . "  " . $gene_name . "  " . $product_name);
            print OUTFILE ("  " . $species . "\n");
            $aa_read_toggle = 1;
        }
        else 
        {
            print OUTFILE ($gi_number .  "   Acc. no.: " . $acc_number . "  " . $gene_name . "  " . $product_name . "  ");
            print OUTFILE ("gi|" . $gi_number . "  " . $species . "\n");
            $aa_read_toggle = 1;
        }
    }
    elsif ($_ =~ /[a-zA-Z]+/) 
    {
        if ($aa_read_toggle == 1) 
        {
            $orf_seq_line = $_;
            chomp($orf_seq_line);
            $orf_seq_line =~ s/[\d]*//g;
            $orf_seq_line =~ s/[\s]*//g;
            $orf_seq_line =~ s/[\W]*//g;
            $orf_seq_line =~ tr/a-z/A-Z/;
            print OUTFILE ($orf_seq_line . "\n");
        }
    }
    elsif ($_ =~ /\\/) 
    {
        if ($aa_read_toggle == 1) 
        {  
            $aa_read_toggle = 0;
        }
    }
}
