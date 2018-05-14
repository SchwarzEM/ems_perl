#!/usr/bin/perl -w

# Program: genpept_to_blast_scripts.pl
#    Erich Schwarz, 8/30/99
# 
# Purpose: Once you've got a bunch of homologs in a folder,
#    how to do batch BlastP or psi-Blast with them?
#    With free-standing blastall and blastpgp from NCBI
#    and the shell scripts produced by this script.

# 1. If not given a genpept file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input GenPept file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$out_blastp_script = ($infile . ".blastp.sh");
$out_psi_blast_script = ($infile . ".psi-blast.sh");
print "The input file is $infile \n";
print "    the output BlastP script is $out_blastp_script \n";
print "    and the output psi-Blast script is $out_psi_blast_script \n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "GenPept file $infile not found\n";
open (OUTBLASTPSCRIPT, ">$out_blastp_script") 
    || die "Couldn't open file $out_blastp_script\n";
open (OUTPSIBLASTSCRIPT, ">$out_psi_blast_script") 
    || die "Couldn't open file $out_psi_blast_script\n";

# 3. Give each script a header that will make it work on my machine
#       and even be intelligible.

$date = `date`;

print OUTBLASTPSCRIPT ('#!/bin/bash '."\n");
print OUTBLASTPSCRIPT ('# '."$out_blastp_script -- script for ");
print OUTBLASTPSCRIPT ("batch BlastP on-site search \n");
print OUTBLASTPSCRIPT ('# Extracted from '."$infile on $date \n");

print OUTPSIBLASTSCRIPT ('#!/bin/bash '."\n");                                                        
print OUTPSIBLASTSCRIPT ('# '."$out_psi_blast_script -- script ");
print OUTPSIBLASTSCRIPT ("for batch psi-Blast on-site search \n");
print OUTPSIBLASTSCRIPT ('# Extracted from '."$infile on $date \n");

# 4. Extract everything I want in a format which is of real use.

@gi_nos = ();

while (<INFILE>) {
	if ($_ =~ /^PID[^0-9]*([0-9]+)$/) {
		$gi_number = $_;
		chomp($gi_number);
		$gi_number =~ s/[^0-9]//g;
		push (@gi_nos , $gi_number);
		}
}

@gi_nos = sort(@gi_nos);

$i = 0;

for ($i = 0; $i <= $#gi_nos; $i++) {	

print OUTBLASTPSCRIPT ("blastall -p blastp -d nr -i $gi_nos[$i] -o $gi_nos[$i]");
print OUTBLASTPSCRIPT (".blastp -e 0.05 -b 0 -I T ;\n");

print OUTPSIBLASTSCRIPT ("blastpgp -d nr -i $gi_nos[$i] -o $gi_nos[$i]");
print OUTPSIBLASTSCRIPT (".psi-blast -e 0.001 -j 3 -I T -F T -h 0.001 -b 0 ;\n");

# e.g.:
# print OUTFILE ("blastall -p blastp -d nr -i 4545123 -o 4545123.blastp -e 0.05");

}

system "chmod a+x $out_blastp_script";
system "chmod a+x $out_psi_blast_script";
