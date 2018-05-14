#!/usr/bin/perl -w

# Program: nrdb_to_clean_tfa.pl
#    Erich Schwarz, 3/20/00
# 
# Purpose: Holm's nrdb.pl script has a few problems with its output:
#             1. It needs "\n" to replace "\\\\*\n".
#             2. It could stand to do "\n" after 68 characters of a line lacking ">".

if ($#ARGV != 0) {
	print "Required: input nrdb file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".cleaned");
print "The input file is $infile; the output file $outfile\n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "GenPept file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Remove crud, break up lines.

while (<INFILE>)
{
    if ($_ =~ />/)
    {
        $_ =~ tr/[\\]{2}[.]+\n/\n/;
        print OUTFILE $_;
    }
    else 
    {
        print OUTFILE $_;

#        $line_to_break = $_;
#            while ($line_to_break =~ /^([\S]{68})/) 
#            {
#                print OUTFILE $1;
#                $line_to_break =~ tr/^([\S]{68})//;
#            }

    }
}
