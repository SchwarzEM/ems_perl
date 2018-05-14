#!/usr/bin/perl -w 

# script_spec_blastps.pl
# Erich Schwarz, 5/9/02

# generate, from a list of sequences, script commands with this format:

# nice ../blastpgp -d ../ATH1.pep -F T -a 2 -I T -e 1e-02 -j 1 -v 1000 -b 1000 -i 9961345_aa_5-370 -m 6 -o 9961345_aa_5-370.psi-blast-x1.ATH1.pep ;
# nice ../blastpgp -d ../Homo.genscan -F T -a 2 -I T -e 1e-02 -j 1 -v 1000 -b 1000 -i 9961345_aa_5-370 -m 6 -o 9961345_aa_5-370.psi-blast-x1.Homo.genscan ;
# nice ../blastpgp -d ../Mus.genscan -F T -a 2 -I T -e 1e-02 -j 1 -v 1000 -b 1000 -i 9961345_aa_5-370 -m 6 -o 9961345_aa_5-370.psi-blast-x1.Mus.genscan ; 

print "List of sequences: ";
chomp ($input = <STDIN>);
$output = $input . ".spec_blast.sh";

@database = qw( ATH1.pep Homo.genscan Mus.genscan fugu_proteome );

open (SEQ_LIST, "$input") || die;
open (BLAST_SHELL, ">$output") || die;

while (<SEQ_LIST>) 
{
    chomp ($input_sequence = $_);
    $i = 0;
    while ($i < 4) 
    {
        print BLAST_SHELL 'nice ../blastpgp -d ../';
        print BLAST_SHELL "$database[$i]";
        print BLAST_SHELL " ";
        print BLAST_SHELL '-F T -a 2 -I T -e 1e-02 -j 1 -v 1000 -b 1000 -i ';
        print BLAST_SHELL "$input_sequence";
        print BLAST_SHELL " -m 6";
        print BLAST_SHELL " -o ";
        print BLAST_SHELL "$input_sequence";
        print BLAST_SHELL ".psi-blast-x1.";
        print BLAST_SHELL "$database[$i]";
        print BLAST_SHELL " ;\n";
        $i += 1;
    }
}

close SEQ_LIST;
close BLAST_SHELL;

system "chmod +x $output";
