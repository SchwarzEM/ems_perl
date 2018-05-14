#!/usr/bin/perl

# sort_annot_generefs.pl; Erich Schwarz, <emsch@its.caltech.edu>, 3/24/04.
# Thanks to Juancarlos Chan for donated Xref code.
# Purpose: convert human-readable annots. into sorted nonredundant normalized gene/reference tables
# Usage: sort_annot_generefs.pl <[input_file(s) or <STDIN>] >[output_files or <STDOUT>] ; can be used in pipeline.
# Input example 1: Provisional_description "rrf-1" Paper_evidence "[cgc4137]" [or] Paper_evidence "[pmid12478294]"
# Input example 2: Provisional_description "T02D1.5" Paper_evidence "[cgc4103]" [or] PMID_evidence "12478294"

use strict;
use LWP::Simple;

my $input_file_path = "";
my $filename = "";
my $gene_label = "";
my %cgcHash;    # hash of cgcs, values pmids
my %pmHash;     # hash of pmids, values cgcs
my %other2loc_hash;  # hash from    other names to locus names
my %seq2loc_hash;    # hash from sequence names to locus names

&populateXref();
&list_seq2loc();
&list_other2loc();

open SORT_UNIQ, "| sort | uniq" or die "Can't pipe to \"|sort | uniq\": $!\n";

foreach $input_file_path (@ARGV) {
    my %genes = "";
    open INPUT_FILE, "$input_file_path" or die "Can't open $input_file_path: $!\n";
    $filename = $input_file_path;
    if ($input_file_path =~ /.*\/+[^\/]+$/) {
        $filename =~ s/.*\/+([^\/]+)$/$1/;
    }
    my @stat_value = stat($input_file_path);
    my @ldate = localtime $stat_value[9];
    foreach my $input_file_line (<INPUT_FILE>) {
        if ($input_file_line =~ /^\w+_description\s+\"([^\"]+)\"\s+Paper_evidence\s+\"\W+(\w+)\W+\"/) {
            $gene_label = $1;
            if (&fix_alt_name($gene_label)) {
                $gene_label = &fix_alt_name($gene_label);
            }
            my $wb_ref  = $2;
            if ($wb_ref =~ /^pmid/ && &checkNumber($wb_ref)) {
                $wb_ref = "cgc" . &checkNumber($wb_ref);
            }
            &output_stats($gene_label, $wb_ref, $ldate[5], $ldate[4], $ldate[3], $filename);
        }
        elsif ($input_file_line =~ /^\w+_description\s+\"([^\"]+)\"\s+PMID_evidence\s+\"(\d+)\"/) {
            $gene_label = $1;
            if (&fix_alt_name($gene_label)) {
                $gene_label = &fix_alt_name($gene_label);
            }
            my $wb_ref  = "pmid" . $2;
            if ($wb_ref =~ /^pmid/ && &checkNumber($wb_ref)) {
                $wb_ref = "cgc" . &checkNumber($wb_ref);
            }
            &output_stats($gene_label, $wb_ref, $ldate[5], $ldate[4], $ldate[3], $filename);
        }
    }
    close INPUT_FILE;
}

close SORT_UNIQ;

sub output_stats {
    # $gene_label, $wb_ref, $ldate[5], $ldate[4], $ldate[3], $filename
    my ($glabel, $ref, $year, $month, $day, $name) = @_;
    if ($glabel =~ /.+[a-zA-Z]{1}/) {
        # trim suffix a-z/A-Z from names
        $glabel =~ s/(.+)[a-z]{1}$/$1/; 
        # deleting A-Z suffix munges eif-3.B->K!
    }
    print  SORT_UNIQ ($year+1900);
    printf SORT_UNIQ "%02d", $month+1;
    printf SORT_UNIQ "%02d", $day;
    print  SORT_UNIQ "\t$name\t$glabel\t$ref\t";
    # name -> etc. order: groups same-dated files in output -- important!
    printf SORT_UNIQ "%02d", $month+1;
    print  SORT_UNIQ "\/";
    printf SORT_UNIQ "%02d", $day;
    print  SORT_UNIQ "\/";
    print  SORT_UNIQ ($year+1900);
    print  SORT_UNIQ "\n";
}

sub list_other2loc {
    # Note!  This is a total kludge, depending on hand-made Tablemaker output; replace with aceserver request!
    open LOC2OTHER, "locus2othername_WS121.txt" or die "Can't open locus2othername_WS121.txt\n";
    foreach my $loc2other_line (<LOC2OTHER>) {
        chomp $loc2other_line;
        $loc2other_line =~ s/\"//g;
        if ($loc2other_line =~ /(\S+)\t(\S+)\s*/) {
              $other2loc_hash{$2} = $1;
        }
    }
}

sub list_seq2loc {
    # Note!  This is a total kludge, depending on hand-made Tablemaker output; replace with aceserver request!
    open LOC2SEQ, "locus2sequence_WS121.txt" or die "Can't open locus2sequence_WS121.txt\n";
    foreach my $loc2seq_line (<LOC2SEQ>) {
        chomp $loc2seq_line;
        $loc2seq_line =~ s/\"//g;
        if ($loc2seq_line =~ /(\S+)\t(\S+\d)[a-zA-Z]/) {
              $seq2loc_hash{$2} = $1;
        }
        elsif ($loc2seq_line =~ /(\S+)\t(\S+)/) {
              $seq2loc_hash{$2} = $1;
        }
    }
}

sub fix_alt_name {
    my $alt_name = shift;
    if ($other2loc_hash{$alt_name}) {
        $other2loc_hash{$alt_name};
    }
    elsif ($seq2loc_hash{$alt_name}) {
        $seq2loc_hash{$alt_name};
    }
}

sub populateXref {
    my $page = get "http://minerva.caltech.edu/~postgres/cgi-bin/cgc_pmid_xref.cgi";
    my @lines = split/\n/, $page;
    foreach my $line (@lines) {
        $line =~ m/<TR><TD ALIGN=CENTER>cgc(\d+)<\/TD><TD ALIGN=CENTER>pmid(\d+)<\/TD><\/TR>/;
        $cgcHash{$1} = $2;
        $pmHash{$2} = $1;
    }
}

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
}

# from the table, generate this second table:
#    date  filename  gene  num_refs  qual_refs
#
# from the table, generate these numbers:
#    date  filename  num_genes  num_refs  mean_refs_per_gene  num_genes_qual_1.0
