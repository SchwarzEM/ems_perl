#!/usr/bin/perl -w

# Program: blastpgp_format-6_to_tfas.pl
#    Erich Schwarz, 5/30/01
# 
# Purpose: 
# 
# psi-BLAST produces hits in protein databases.  Often it converges 
# on a broad general consensus from many homologs, that does not fully 
# capture the diversity of the more far-flung members of a protein 
# family; redoing psi-BLAST with individual hits in a "converged" set 
# can detect yet more homologs at good stringency that were missed 
# by the overall consensus.
# 
# Proteins are modular.  In a large protein family, it is entirely possible 
# for a divergent homolog to exist whose family-homologous domain is fused 
# to another domain with entirely different function and homologies.
# 
# A simple solution to psi-BLAST's general consensus -- iterating 
# psi-BLAST with individual hits in a final set -- will founder if 
# even a few of these individual hits are chimeric; there will arise 
# large numbers of spurious homologs, and sorting out which are which 
# by hand defeats the whole point of automatic searching.
# 
# A straightforward solution to this is to extract only those regions of 
# the final homolog set, from the psi-BLAST report, that are actually 
# detected by the first query sequence (and thus have real homology to it).
# Using these extracted sequences as seeds for further psi-BLAST searches 
# may lower sensitivity, but is at least bound to prevent chimeric proteins 
# from sabotaging iterative psi-BLAST.  Better yet, any new hits that are 
# detected in this manner should themselves be extractable as limited regions
# of homology from the second set of psi-BLAST reports.
# 
# This script is designed to work with psi-BLAST reports generated on Linux 
# by the freestanding "blastpgp" binary provided by NCBI in early 2001.  The 
# following arguments need to be used in order to have the report be properly 
# extractable by this script:

# blastpgp ... -I T -e [X^1] -h [X^1] -j [X^2] -v [X^3] -b [X^3] -m 6 ...
# 
# The key thing here is the argument "-m 6".  This produces a specific 
# kind of format in which the homologous regions of hits are specifically 
# shown, residue by residue, which the program needs to generate FASTA 
# files containing the appropriate sequences.
# 
# Other actual values are arbitrary but: -v and -b should be identical and 
# large; -j should be large enough to induce convergence (the program 
# looks for "CONVERGED!" as a sign to start reading the output); -e 
# and -h should be identical and acceptably stringent.  For example:
# 
# blastpgp ... -I T -e 1e-06 -h 1e-06 -j 20 -v 1000 -b 1000 -m 6 ...

# This is the sort of output that the program will work 
# on (very briefly excerpted):

# ...
# Results from round 5
# 
# 
#                                                                    Score     E
# Sequences producing significant alignments:                        (bits)  Value
# Sequences used in model and found again:
# 
# C31A11.5 CE15656    (HINXTON) TR:O45282 protein_id:CAB05687.1     648  0.0
# F09B9.1 CE15758    (HINXTON) TR:Q19239 protein_id:CAA90058.1      643  0.0
# F52F10.4 CE19884    (ST.LOUIS) TR:Q9UAQ9 protein_id:AAC69230.1    628  e-180
# F52F10.3 CE19883    (ST.LOUIS) TR:Q9UAQ8 protein_id:AAC69231.1    628  e-180
# ...
# C49D10.7 CE08848    (ST.LOUIS) TR:O16601 protein_id:AAC71185.1    150  3e-36
# T08H10.4 CE07477    (ST.LOUIS) TR:Q22351 protein_id:AAA97991.1    100  6e-21
# F47F6.3 CE10708    (ST.LOUIS) TR:P91310 protein_id:AAC71112.1      98  1e-20
# Y69A2AR.11 CE27510    (ST.LOUIS)                                   91  2e-18
# 
# Sequences not found previously or not previously below threshold:
# 
# 
# CONVERGED!
# QUERY      125 YLYNIGKRFLAITQPFDAARN-PDASTLCRRQMHQF---LNALDNF-------------- 166
# C31A11.5   65  ------------------------IEQKCTKEELKT---IE--DNF-------------- 81
# F09B9.1    34  ------------------------LSPQCLNDTQTW---IKSLELFAGVSEACLKKQTCS 66
# F52F10.4   34  -----------------------NLSAQCLNDTNTW---INSLEIFGTLYAECIIMKKCN 67
# F52F10.3   16  -----------------------DLSEKCLNDTDTW---LKSLEIFSTVSVECLTLNNCT 49
# ...
# C49D10.7   28  ---------------------F-FPD--TFP----------NGYLGI---LRSSKFLMTM 50
# T08H10.4   22  ---------------------F-FPK--TFP----------NGYIGVDMFFVLSGFLMAM 47
# F47F6.3        ------------------------------------------------------------
# Y69A2AR.11 31  ---------------------F-YPD--TFP----------NGYLGVDQFFVLSGFLMCM 56
# ...
# QUERY      770 ---------APFQNLK--------KLLIKRP 783
# C31A11.5   652 ---------IPFLKLE--------KMLIE-- 663
# F09B9.1    651 ---------IPTLKLE--------KMLIE-- 662
# F52F10.4   639 ---------VPILKLE--------KMLIEK- 651
# F52F10.3   633 ---------MPILKLE--------KMLIES- 645
# ...
# C49D10.7   139 ------------------------------- 140
# T08H10.4   75  ------------------------------- 76
# F47F6.3    200 ---------KQYLQLD--------FK----- 208
# Y69A2AR.11 65  ------------------------------- 66
#   Database: wormpep51
#     Posted date:  May 21, 2001 12:02 AM
#   Number of letters in database: 8,676,709
#   Number of sequences in database:  19,781
# ...

# Though in an output done on nr, it'll be likely to have hits like this:

# QUERY    358 --PWFT---AF-SLDKNLR------WLF---ST-S--------------S-A--P-G--- 380
# 10728260     ------------------------------------------------------------
# 5815231      ------------------------------------------------------------
# 7291079  260 --PLLI---AF-SVLTNAP------KIF---TV-K--------------K-VNNP-N--- 284
# 7293604  254 --PVVE---AF-SARANSR------ALFRIVDT-K--------------A-N--P-N--- 279

# These are from files headed ">gi|X|...", but apparently blastpgp cleans that 
# up when its outputs are in "-m 6" mode.  Good thing.

# OK, on with the coding.  Basic strategy:

# 1. Open the psi-BLAST report to be read.
# 2. Scan for "CONVERGED!" in each input line.  Until it is found do nothing.
# 3. Once it's found, become ready to create output arrays.
#
# 4. With each new line after that, 
#        4A. if there are 1+ non-blank initial characters <word>,
#
#                create a list of <name>s, a header array and 
#                  a residue array.  Each array should be named <word>_header 
#                  or <word>_residue.  The <name>s list should eventually give all
#                  names found.
#                The description array will hold <word>, <lowest residue>,
#                  <highest residue>, and <database>.  Start with values "0".
#                The residue array will only have one entry, <word>, at first.
#
#        4B. else, stop creating output arrays and get ready to fill them.
#
# 5. With each new line after that,
#        5A. If "Database:" is in the input line:
#
#                take the database name
#                insert it into slot 4 of all header arrays
#                  by reading them off the name list
#                go to printout of array contents
#
#        5B. If a word starts the line:
#
#                open the header array named that word
#                open the residue array named that word, and:
#                  take the leftmost residue number, if there, and insert it into the 
#                    starting-residue slot, iff the value isn't 0 AND it's lower
#                    than the value already there.
#                  take the rightmost residue number, if there, and insert it into 
#                    the ending-residue slot, iff it's higher than the present value.
#                spit all the amino-acids in the rest of the line after that
#                  into the residue array.
# 
#                close the header and residue arrays named that word
# 
# 6. To print out array contents:
#
#        6A. Create a *.tfa file named after <name>
#        6B. Print ">" followed by <name> [header array $1]
#                  "residues <start>-<end> -- database: <database> \n"
#        6C. Print residues from residue array.  Try to print residues N*60 + 1-60,
#                     where N=0 
#                 if it fails, print "\n" and stop.
#                 if it succeeds, increment N by +1 and repeat attempt.
#        Eventually this should give a 60-residue/line FASTA file with name, residues,
#                 and database.
#
# 7. To print out shell script for doing iterative psi-Blast:
#        read each <name> in turn into a pre-written line command, written to 
#        simply allow database querying (j=0, b=0) from each new sequence.

# Beginning of code:

# 1. Open the psi-BLAST report to be read.

if ($#ARGV != 0) {
	print "Required: input psi-BLAST output.\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".tfa");
print "The input file is $infile.\n";

# Initialize control variables to prevent generation of arrays 
# or reading of sequences before the converged set.

$array_write_toggle = 0;
$sequence_read_toggle = 0;

open (INFILE, "$infile") || die "psi-BLAST output $infile not opened\n";

# 2. Scan for "CONVERGED!" in each input line.  Until it is found do nothing.

while (<INFILE>) {
	if ($_ =~ /CONVERGED!/) {
		$array_write_toggle = 1;
		}
	elsif ($_ =~ /^PID[^0-9]*([0-9]+)$/) {
		$aa_read_toggle = 0;
		$gi_number = $_;
		chomp($gi_number);
		$gi_number =~ s/[^0-9]//g;
		$species = "no_species_name";
		$gene_name = "no_gene_name";
		$product_name = "no_product_name";
		print OUTFILE (">");
	}
	elsif ($_ =~ /ORGANISM[\s]+[.]*/) {
		$species = $_;
		chomp($species);
		$species =~ s/ORGANISM[\s]+//g;
	}
	elsif ($_ =~ /gene=/) {
		$gene_name = $_;
		chomp($gene_name);
		$gene_name =~ s/[\W]*gene=//;
		$gene_name =~ s/\"//g;
	}
	elsif ($_ =~ /product=/) {
		$product_name = $_;
                chomp($product_name);  
                $product_name =~ s/[\W]*product=//;
                $product_name =~ s/\"//g;
	}

# Direct next output by whether $format_toggle = 1 and whether gene has name.

	elsif ($_ =~ /ORIGIN/) {

	if ($format_toggle == 1) {

	if ($gene_name =~ /no_gene_name/) {
		print OUTFILE ($gi_number . "  " . $gene_name . "  " . $product_name . "  ");
		print OUTFILE ("gi|" . $gi_number . "  " . $species . "\n");
		$aa_read_toggle = 1;
	}
        else {   # I.e., if $gene_name is non-trivial...
                print OUTFILE ($gene_name . "  " . $product_name . "  ");                    
                print OUTFILE ("gi|" . $gi_number . "  " . $species . "\n");
                $aa_read_toggle = 1;
        }
	}

        else {
                print OUTFILE ($gi_number . "  " . $gene_name . "  " . $product_name . "  ");
                print OUTFILE ("gi|" . $gi_number . "  " . $species . "\n");
                $aa_read_toggle = 1;
        }
	}

	elsif ($_ =~ /[a-zA-Z]+/) {
		if ($aa_read_toggle == 1) {
			$orf_seq_line = $_;
			chomp($orf_seq_line);
			$orf_seq_line =~ s/[\d]*//g;
                	$orf_seq_line =~ s/[\s]*//g;
                	$orf_seq_line =~ s/[\W]*//g;
			$orf_seq_line =~ tr/a-z/A-Z/;
			print OUTFILE ($orf_seq_line . "\n");
			}
	}
        elsif ($_ =~ /\\/) {
                if ($aa_read_toggle == 1) {  
                        $aa_read_toggle = 0;
                        }
        }

}
