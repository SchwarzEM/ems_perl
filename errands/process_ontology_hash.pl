#!/usr/bin/perl -w

# Program: process_ontology_hash.pl
#    Erich Schwarz, 6/29/01
# 
# Purpose: Create hash of "GO:---" to "P" for public gene 
# assignment table.  We want just one line per GO term.

print "What will input process ontology file be?\n";
print "Default is \"process.ontology\": ";
$infile = <STDIN>;
chomp ($infile);
unless ($infile =~ /\w+/)
{
    $infile = "process.ontology";
}

$outfile = ($infile . ".hash");
print "The input file is $infile; the final output file $outfile\n";

# 2. Open the files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found\n";
open (OUTFILE, ">$outfile.1") || die "Couldn't open file $outname.1\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
    $line_to_be_extracted = $_;
    {
        while ($line_to_be_extracted =~ /(GO:\d{7})/) 
        {
            print OUTFILE "$1\n";
            $line_to_be_extracted =~ s/(GO:\d{7})//;
	}
    }
}

close INFILE;
close OUTFILE;

# Output file needs to be sort-ed and uniq-ed; it'd be 
#    easiest to do this by system calls.

system "sort $outfile.1 > $outfile.2";
system "uniq $outfile.2 > $outfile.3";

open (INFILE, "$outfile.3") || die "Wormpep file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname.1\n";

while (<INFILE>)  
{
    print OUTFILE $_;
    print OUTFILE "P\n";
}

# Cleanup:

system "rm $outfile.1; rm $outfile.2; rm $outfile.3";

# Note that to get this back in as a hash it has to be 
#   treated as a list.
