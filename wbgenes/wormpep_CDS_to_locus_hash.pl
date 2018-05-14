#!/usr/bin/perl

# Program: wormpep_CDS_to_locus_hash.pl
#    Erich Schwarz, 1/24/02.
# 
# Purpose: Get CDS-to-locus hash out of wormpep headers.
#   (This works for wormpep53 and seems to still be working for wormpep72;
#    note, though, that it might break if the wormpepN header format 
#    seriously changed.)

# 1. If not given a wormpep file as argument, ask for its name.

if ($#ARGV != 0) 
{
    print "What will input wormpepN file be? ";
    $infile = <STDIN>;
} else 
{
    $infile = $ARGV[0];
}
chomp ($infile);
$outfile = ("wormpep.CDS-to-locus.hash");
print "The input file is $infile; the output file $outfile\n";

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname\n";

# 3. Extract everything I want in a format which is of real use.

# Note: "." is changed to "_" in CDSes to get around 
#   the apparent choking of hashes by "." in scalars.
#   This could probably be handled in some more clever way, 
#   but not at my current Perl skill level just yet.

while (<INFILE>) 
{
    if ($_ =~ /^>(\w+\.\w+) .+ locus:(\S+)\s+/) 
    {
        $cds = $1;
        $locus = $2;
        $cds_hash_friendly = $cds;
        $cds_hash_friendly =~ s/\./\_/g;
        print OUTFILE "$cds_hash_friendly\n";
        print OUTFILE "$locus\n";
    }
    elsif ($_ =~ /^>(\w+\.\w+) /)
    {
        $cds = $1;
        $cds_hash_friendly = $cds;
        $cds_hash_friendly =~ s/\./\_/g;
        print OUTFILE "$cds_hash_friendly\n";
        print OUTFILE "$cds\n";
    }
    elsif ($_ =~ /^>/)
    {
        print "WARNING: The input line \"$_\" could not be interpreted by this script.\n";
    }
}
