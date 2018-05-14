#!/usr/bin/perl -w 

# script_spec_TBlastN.pl
# Erich Schwarz <emsch@its.caltech.edu>, 12/04/02

# Purpose: generate, from lists of sequences and protein databases, scripted TBlastNs for each versus each.

$input = "";
$list_file = "";
$threshold = "";
@database = "";
$indiv_database = "";

print "List of sequences: ";
chomp ($input = <STDIN>);

print "List of protein sets (e.g. -- nr; wormpep91; HumanPep4): ";
chomp ($list_file = <STDIN>);

print "E value to be used as a threshold in the searches: ";
chomp ($threshold = <STDIN>);

open (LIST_FILE, "$list_file");
@database = (<LIST_FILE>);
chomp(@database);

foreach $indiv_database (@database) 
{
    $output = $input . "." . $indiv_database ."_psi-blast-x01.sh";

    open (SEQ_LIST, "$input") || die;
    open (BLAST_SHELL, ">$output") || die;

    while (<SEQ_LIST>) 
    {
        chomp ($input_sequence = $_);
        {
            # Desired format: something like "DRD.psi-blast-x01.nr_1e-06"

            print BLAST_SHELL 'nice ../blastall -p tblastn -d ../';
            print BLAST_SHELL "$indiv_database ";
            print BLAST_SHELL "-F T -a 2 -I T ";
            print BLAST_SHELL "-e ";
            print BLAST_SHELL "$threshold ";
            print BLAST_SHELL "-v 1000 -b 1000 -i ";
            print BLAST_SHELL "$input_sequence";
            print BLAST_SHELL " -m 6";
            print BLAST_SHELL " -o ";
            print BLAST_SHELL "$input_sequence";
            print BLAST_SHELL ".TBlastN.";
            print BLAST_SHELL "$indiv_database";
            print BLAST_SHELL "_";
            print BLAST_SHELL "$threshold";
            print BLAST_SHELL " ;\n";
        }
    }
    close SEQ_LIST;   
    close BLAST_SHELL;
    system "chmod +x $output";
}
