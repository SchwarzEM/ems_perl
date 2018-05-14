#!/usr/bin/perl

# acc2gpep.pl: Erich Schwarz <emsch@its.caltech.edu>, 3/7/06; revised 5/12/08.
# Purpose: given accession numbers, get GenPept text while complying with NCBI guidelines.
# READ: http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html#UserSystemRequirements
# Documentation:
#     http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html
#     http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=coursework.section.brief
#     http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html
# Originally derived from "Pham, Vyvy (NIH/NLM/NCBI)" <pham@ncbi.nlm.nih.gov>, 5/12/04.

use strict;
use warnings;
use LWP::Simple;

my $q = q{};
my $i = 0;

while (my $input = <>) {
    $i++; 
    chomp $input; 
    if ( $input !~ / \A \W* (\w+) \W* \z /xms ) { 
        warn "Ambiguous: filter or transform?:  $input\n"; 
    }
    $q .= "$input,";
    if ($i == 100) {
        $i = 0; 
        post();
        sleep 4;  # Pace requests: 1x per >= 3 seconds.
        $q = q{};
    }
}

post();

sub post() { 
    my $fullquery = 'http://eutils.ncbi.nlm.nih.gov/'                 .
                    'entrez/eutils/efetch.fcgi?'                      .
                    'db=protein&retmode=text&rettype=gp&retmax=10000' .
                    '&id=' . "$q" . '&email=emsch@its.caltech.edu';      # Use "email" so NCBI can reply.
    my $d = get("$fullquery");
    print "$d\n";
} 

