#!/usr/bin/perl

# fix_WBaseGO_refs.pl
# Erich Schwarz, <emsch@its.caltech.edu>, 3/22/04.
# Thanks to Juancarlos Chan for donated Xref code (slightly changed here), and to Carol Bastiani for pointing out bug.
# Purpose: convert WB:[cgc1470] or WB:[pmid1782862], both of which are dorky, into WB:[cgc1470]|PMID:1782862 for gene_assocations.wb.
# Usage: fix_WBaseGO_refs.pl <[input_file or <STDIN>] >[output_file or <STDOUT>] ; can be used in pipeline.

use LWP::Simple;

my $input_line = "";

my %cgcHash;    # hash of cgcs, values pmids
my %pmHash;     # hash of pmids, values cgcs
&populateXref();

while (<>) {
    chomp($input_line = $_); 
    if ($input_line =~ /(\t[^\t]*PMID:(\d+)[^\t]*\t)/) {
        my $input_pmid_no = $2;
        my $input_pmid_word = "pmid" . $input_pmid_no;
        if (&checkNumber($input_pmid_word)) {
            print "$`";
            print "\tWB:[cgc";
            print &checkNumber($input_pmid_word);
            print "]|PMID:";
            print $input_pmid_no;
            print "\t$'\n";
        }
        else {
            print "$`";
            print "\tPMID:";
            print $input_pmid_no;
            print "\t$'\n";
        }
    }

    elsif ($input_line =~ /(\t[^\t]*WB:\[(cgc\d+)\][^\t]*\t)/) {
        my $input_cgc = $2;
        if (&checkNumber($input_cgc)) {
            print "$`";
            print "\tWB:[";
            print $input_cgc;
            print "]|PMID:";
            print &checkNumber($input_cgc);
            print "\t$'\n";
        }
        else {
            print "$`";
            print "\tWB:[";
            print $input_cgc;
            print "]\t$'\n";
        }
    }

    elsif ($input_line =~ /(\t[^\t]*WB:\[pmid(\d+)\][^\t]*\t)/) {
        my $input_pmid_no = $2;
        my $input_pmid_word = "pmid" . $input_pmid_no;
        if (&checkNumber($input_pmid_word)) {
            print "$`";
            print "\tWB:[cgc";
            print &checkNumber($input_pmid_word);
            print "]|PMID:";
            print $input_pmid_no;                                                                                                                    
            print "\t$'\n";
        }
        else {
            print "$`";
            print "\tPMID:";
            print $input_pmid_no;
            print "\t$'\n";
        }
    }

    elsif ($input_line =~ /\t[^\t]*WB:\[\S+\][^\t]*\t/) {  # fallback position -- just print the line
        print "$input_line\n";
    }
    else { 
        print "$input_line\n";
        warn "Possible non-referenced gene association\n";
    }
}

sub populateXref {
    # if not found, get ref_xref data to try to find alternate
    my $page = get "http://minerva.caltech.edu/~postgres/cgi-bin/cgc_pmid_xref.cgi";
    my @lines = split/\n/, $page;
    foreach my $line (@lines) {
        $line =~ m/<TR><TD ALIGN=CENTER>cgc(\d+)<\/TD><TD ALIGN=CENTER>pmid(\d+)<\/TD><\/TR>/;
        $cgcHash{$1} = $2;
        $pmHash{$2} = $1;
    }
} # sub populateXref

sub checkNumber {
    my $number_coming_in = shift;
    # if -- a capitalization insensitive pmid with possibly brackets around it
    if ($number_coming_in =~ m/\[?[pP][mM][iI][dD](\d+)\]?/) {
        if ($pmHash{$1}) {
            $pmHash{$1};          # returns this value -- no 'print' needed!
        }
    }
    # elsif -- a capitalization insensitive cgc with possibly brackets around it
    elsif ($number_coming_in =~ m/\[?[cC][gG][cC](\d+)\]?/) {
        if ($cgcHash{$1}) {
            $cgcHash{$1};         # again, 'print' is worse than useless here
        }
    }
    # else -- neither a cgc nor pmid, which shouldn't be fed to this subroutine!
    else { 
        die "Useless invocation of populateXref subroutine!\n";
    }
} # sub checkNumber

