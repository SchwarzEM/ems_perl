#!/usr/bin/perl -w 

# mutate_text.pl
# Erich Schwarz, 11/21/01.
#
# Very simple minded: convert strings in a text file
# into different strings.  Of use in converting Blast
# scripts quickly into different ones.

print "Script to mutate?:  ";
$infile=<STDIN>;
chomp($infile);
print "\n";
print "Name of mutant product script?:  ";
$outfile=<STDIN>;
chomp($outfile);
print "\n";
print "Exact string in first script to alter?:  ";
$instring=<STDIN>;
chomp($instring);
print "EXACT altered string to use in second text?:  ";
$outstring=<STDIN>;
chomp ($outstring);

open INFILE, "$infile" || die "Can't open $infile! $!\n";
open OUTFILE, ">$outfile" || die "Can't open $outfile! $!\n";

while (<INFILE>) 
{
    $inline = $_;
    $inline =~ s/$instring/$outstring/g;
    print OUTFILE $inline;
}

close INFILE;
close OUTFILE;
