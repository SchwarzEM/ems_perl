#!/usr/bin/perl -w

# Program: tfa_split.pl
# Erich Schwarz <emsch@its.caltech.edu>, 9/25/03
# 
# Purpose: split .tfa files with multiple entries into many files with single entries; make blastp and extractseq scripts for them.

# Note: adapted for dual-processor Pentium! Delete argument "-a 2" for single-processors.

# 1. Open a .tfa file

if ($#ARGV != 0) 
{
    print "Required: input multiple FASTA file!\n";
    print "What will input file be? ";
    $infile = <STDIN>;
} 
else 
{
    $infile = $ARGV[0];
}
chomp ($infile);
$outfolder = ($infile . "_indiv_tfas_folder");

print "The input file is $infile; the output *directory* $outfolder\n";

# Open infile and create outfolder and array for script later.

open (INFILE, "$infile") || die "Mass-FASTA file $infile not found\n";
mkdir ($outfolder,0777) || die "Couldn't create directory $outfolder\n";
chdir ($outfolder) || die "Can't change to directory $outfolder\n";
@protein_names = ();

# 2. Scan to ">", then spit out a series of appropriately named files.

while (<INFILE>) 
{
    if ($_ =~ />(\S+)\s+/) 
    {
        close (OUTFILE);
        open (OUTFILE, ">$1") || die "Can't open outfile\n"; #line40
        push (@protein_names , $1);
        print OUTFILE $_;
    } 
    else 
    {
        print OUTFILE $_; 
    }
}

# 3. and return to original directory.

close (OUTFILE);

# chdir ("..") || die "Can't return to original directory\n";

$outscript  = ($infile . "_blastp-script");
$outscript2 = ($infile . "_extractseq-script");
$outfile2   = ($infile . "_5000bp");

open (OUTFILE, ">$outscript") || die "Can't open outfile\n";
open (OUTFILE2, ">$outscript2") || die "Can't open outfile\n";

@protein_names = sort(@protein_names);
$i = 0;
for ($i = 0; $i <= $#protein_names; $i++) 
{
    print OUTFILE ("nice ../blastpgp -d ../nr -i $protein_names[$i] -o $protein_names[$i]");
    print OUTFILE (".gap-blast -e 0.001 -F T -b 0 -I T -a 2 ;\n");
}

# Also want to be able to run trimming with extractseq from EMBOSS.
# cat glh-1 | extractseq -filter -regions "1-200" >> test

$i = 0;
for ($i = 0; $i <= $#protein_names; $i++)
{
    print OUTFILE2 ("nice cat $protein_names[$i] | extractseq -filter -regions \"1-5000\" >> $outfile2 ;\n");
}

# This would be what I would use for ungapped BlastP, given past experience:
#    print OUTFILE ("blastall -p blastp -d nr -i $protein_names[$i] -o $protein_names[$i]"); #line60
#    print OUTFILE (".blastp -e 0.001 -b 0 -I T ;\n");

# In some cases the shell doesn't know what to do with a local command: insert "./ ..."
# But I want to run this inside a directory *below* blastpgp: insert "../ [...] ../nr [...]"

close (OUTFILE);

system "chmod a+x $outscript";
system "chmod a+x $outscript2";

$outlist = ($infile . "_gi-list");
open (OUTFILE, ">$outlist") || die "Can't open outfile\n";
$i = 0;

for ($i = 0; $i <= $#protein_names; $i++)
{
    print OUTFILE ("$protein_names[$i]\n");
}

close (OUTFILE);

chdir ("..") || die "Can't return to original directory\n";

