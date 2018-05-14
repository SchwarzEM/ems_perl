#!/usr/bin/perl

# update_KOG_gi_nos.pl
# from "Pham, Vyvy (NIH/NLM/NCBI)" <pham@ncbi.nlm.nih.gov>, 5/12/04
# intended substrate: ftp://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kyva=gb

use LWP::Simple;

while (<>) {
    $i++; 
    chomp; 
    split;
    $q .= "$_[1],";
    if ($i == 500) {
        $i = 0; 
        post(); 
        $q = "";
    }
}

post();

sub post() { 
    $d = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&retmode=xml&retmax=1000&id=$q");
    print "$d\n";
} 
