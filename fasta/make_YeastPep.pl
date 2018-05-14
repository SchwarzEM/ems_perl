#!/usr/bin/perl -w

# make_YeastPep.pl -- v. 2, 1/2/01
# 
# Purpose: make "YeastPep" from yeast-orf_trans.fasta
# ...............................................
# Input header examples:
# 
# >ORFP:Q0010 Q0010, Chr Mito from 3952-4338
# >ORFP:YAL001C TFC3, Chr I from 147591-147660,147751-151163, reverse complement
# >ORFP:YAL009W SPO7, Chr I from 135852-136631
#
# Output header examples:
# 
# >Q0010 Q0010
# >YAL001C TFC3
# >YAL009W SPO7
# ...............................................

# 1. Open yeast-orf_trans.fasta or its equivalent.

if ($#ARGV != 0) {
        print "Required: yeast-orf_trans.fasta or its equivalent.\n";
        print "What will input file be? ";
        $infile = <STDIN>;
} else {
        $infile = $ARGV[0];
}
chomp ($infile);
$outfile = ("YeastPep_" . $infile);

print "The input file is $infile; the output file $outfile\n";

open (INFILE, "$infile") || die "Mass-FASTA file $infile not found. $! \n";
open (OUTFILE, ">$outfile") || die "Can't open outfile $outfile. $! \n";

# 2. Type protein normally, but edit header (">") lines, 
#    and expunge "*" stop residues from one-letter protein sequence.

while (<INFILE>) 
{
    if ($_ =~ /^>ORFP:(\w+)\s+(\w+)/) 
    {
        print OUTFILE (">" . $1 . " " . $2 . "\n");
    }
    elsif ($_ =~ /^(\w+)\*/)
    {
        print OUTFILE ($1 . "\n");
    }
    else 
    {
        print OUTFILE ($_);
    }
}
