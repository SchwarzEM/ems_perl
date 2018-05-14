#!/usr/bin/perl -w 

# script_spec_tblastns.pl
# Erich Schwarz, 5/9/02

# generate, from a list of sequences, script commands with this format:

# ../est_human.Z  ../est_mouse.Z  ../est_others.Z

# nice ../blastall -p tblastn -d ../est_human  -F T -a 2 -I T -e 1e-02 -v 1000 -b 1000 -i 18588797_aa_2-290 -m 6 -o 18588797_aa_2-290.tblastn.est_human ;
# nice ../blastall -p tblastn -d ../est_mouse  -F T -a 2 -I T -e 1e-02 -v 1000 -b 1000 -i 18588797_aa_2-290 -m 6 -o 18588797_aa_2-290.tblastn.est_mouse ;
# nice ../blastall -p tblastn -d ../est_others -F T -a 2 -I T -e 1e-02 -v 1000 -b 1000 -i 18588797_aa_2-290 -m 6 -o 18588797_aa_2-290.tblastn.est_others ;

print "List of sequences: ";
chomp ($input = <STDIN>);
$output = $input . ".spec_tblastn.sh";

@database = qw( est_human est_mouse est_others );

open (SEQ_LIST, "$input") || die;
open (BLAST_SHELL, ">$output") || die;

while (<SEQ_LIST>) 
{
    chomp ($input_sequence = $_);
    $i = 0;
    while ($i < 3) 
    {
        print BLAST_SHELL 'nice ../blastall -p tblastn -d ../';
        print BLAST_SHELL "$database[$i]";
        print BLAST_SHELL " ";
        print BLAST_SHELL '-F T -a 2 -I T -e 1e-02 -v 1000 -b 1000 -i ';
        print BLAST_SHELL "$input_sequence";
        print BLAST_SHELL " -m 6";
        print BLAST_SHELL " -o ";
        print BLAST_SHELL "$input_sequence";
        print BLAST_SHELL ".tblastn.";
        print BLAST_SHELL "$database[$i]";
        print BLAST_SHELL " ;\n";
        $i += 1;
    }
}

close SEQ_LIST;
close BLAST_SHELL;

system "chmod +x $output";
