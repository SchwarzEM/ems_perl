#!/usr/bin/perl

# genpept_to_clean_tfa.pl: Erich Schwarz <emsch@its.caltech.edu>, 3/8/2006
# Purpose: Get useful FASTA files with informative headers from GenPept downloads; use input and output text streams.
# Recent debug: sporadic outputs from NCBI with DNA instead of protein files need to not be read!
# New feature: NCBI TaxIDs now provided like this -- "Shewanella oneidensis [NCBI_TaxID=70863]"

$aa_read_toggle = 0;
$all_read_toggle = 0;

while (<>) {
    if ($_ =~ /^LOCUS\s+\w+\s+\d+\saa\s/) { 
        $all_read_toggle = 1;
        $aa_read_toggle = 0;
    }
    elsif ($_ =~ /^LOCUS/) {
        $all_read_toggle = 0;
        $aa_read_toggle = 0;
    }
    elsif (($all_read_toggle) and ($_ =~ /^ACCESSION\s+(\S+)\s*/)) {    # This is the invariant identifier used by NCBI.
        chomp($acc_number = $1);
        $aa_read_toggle = 0;
        $species = "no_species_name";
        $ncbi_taxno = "none";
        $gene_name = "no_gene_name";
        $product_name = "no_product_name";
        print (">");
    }
    elsif (($all_read_toggle) and ($_ =~ /^VERSION.+GI.([0-9]+)$/)) {    # Format for GI numbers: VERSION     AAF57416.1  GI:7302326
        chomp ($gi_number = $1);
        $aa_read_toggle = 0;
        $gi_number =~ s/[^0-9]//g;
    }
    elsif (($all_read_toggle) and ($_ =~ /ORGANISM[\s]+[.]*/)) {
        chomp($species = $_);
        $aa_read_toggle = 0;
        $species =~ s/ORGANISM[\s]+//g;
        $species =~ s/\<\/a>//;         # this is a very specific kludge, designed in
        $species =~ s/\<a href.+\>//;   # 5/15/02 to remove NCBI cruft from headers...
    }
    elsif (($all_read_toggle) and ($_ =~ /^\s+\/db_xref=\"taxon:(\d+)\"/)) { 
        $ncbi_taxno = '[NCBI_TaxID=' . $1 . ']';
        $aa_read_toggle = 0;
    }
    elsif (($all_read_toggle) and ($_ =~ /gene=/)) {
        chomp($gene_name = $_);
        $aa_read_toggle = 0;
        $gene_name =~ s/[\W]*gene=//;
        $gene_name =~ s/\"//g;
    }
    elsif (($all_read_toggle) and ($_ =~ /product=/)) {
        chomp($product_name = $_);
        $aa_read_toggle = 0;
        $product_name =~ s/[\W]*product=//;
        $product_name =~ s/\"//g;
    }   
    elsif (($all_read_toggle) and ($_ =~ /ORIGIN/)) { 
        print $acc_number;
        print "   gi|" . $gi_number;
        print "  " . $gene_name;
        print "  " . $product_name;
        print "  " . $species;
        print " "  . $ncbi_taxno . "\n";
        $aa_read_toggle = 1;
    }
    elsif (($all_read_toggle) and ($_ =~ /[a-zA-Z]+/)) {
        if ($aa_read_toggle == 1) {
            chomp($orf_seq_line = $_);
            $orf_seq_line =~ s/[\d]*//g;
            $orf_seq_line =~ s/[\s]*//g;
            $orf_seq_line =~ s/[\W]*//g;
            $orf_seq_line =~ tr/a-z/A-Z/;
            print $orf_seq_line;
            print "\n";
        }
    }
    elsif (($all_read_toggle) and ($_ =~ /\/\//)) {
        $aa_read_toggle = 0;
        $all_read_toggle = 0;
    }
}

