#!/usr/bin/perl -w

# create_wormpepN_hashes.pl
#    Erich Schwarz, 1/30/02.
# 
# Purpose: Get CDS-to-locus hash and CDS-to-CEnumber hash
#          out of WormpepN headers.

#   (This works for wormpep53 and seems to still be working for wormpep72;
#    note, though, that it might break if the wormpepN header format 
#    seriously changed.)

# 1. If not given a wormpep file as argument, ask for its name.

print "What will input wormpepN file be? ";
$infile = <STDIN>;

chomp ($infile);

$CDS_to_locus_hash    = "wormpep.CDS-to-locus.hash";
$CDS_to_CEnumber_hash = "wormpep.CDS-to-CEnumber.hash";

print " The input file is                       $infile\n";
print " the output CDS-to-locus hash is         $CDS_to_locus_hash\n";
print " and the output CDS-to-CEnumber hash is  $CDS_to_CEnumber_hash\n";

# 2. Open the WormpepN and output files or die.

open (INFILE, "$infile") || die "Wormpep file $infile not found\n";
open (CDS_TO_LOCUS_HASH, ">$CDS_to_locus_hash") || die "Couldn't open file $CDS_to_locus_hash\n";
open (CDS_TO_CENUMBER_HASH, ">$CDS_to_CEnumber_hash") || die "Couldn't open file $CDS_to_CEnumber_hash\n";

# 3. Extract everything I want in a format which is of real use.

# Note: "." is changed to "_" in CDSes to get around 
#   the apparent choking of hashes by "." in scalars.
#   This could probably be handled in some more clever way, 
#   but not at my current Perl skill level just yet.

while (<INFILE>) 
{

# >AC7.1 CE20434   G-protein coupled receptor status:Predicted TR:Q22876 protein_id:AAB03418.2
# >AC7.2 CE25736  locus:soc-2  status:Predicted TR:Q22875 protein_id:AAB03417.3
# >AC7.3 CE07653    status:Predicted protein_id:AAB03419.1
# >B0350.2C CE06704  locus:unc-44  status:Predicted
# >B0350.2D CE06705  locus:unc-44  status:Predicted
# >B0365.3 CE07721  locus:eat-6 Na(+)\/K(+) ATPase alpha subunit status:Confirmed TR:P90735 protein_id:CAB02694.1



    if ($_ =~ /^>(\S+)\s+(CE\d+)\D{1}.*locus:(\S+)\s+/) 
    {
        $cds = $1;
        $CEnumber = $2;
        $locus = $3;
        $cds_hash_friendly = $cds;
        $cds_hash_friendly =~ s/\./\_/g;
        print CDS_TO_LOCUS_HASH    "$cds_hash_friendly\n";
        print CDS_TO_LOCUS_HASH    "$locus\n";
        print CDS_TO_CENUMBER_HASH "$cds_hash_friendly\n";
        print CDS_TO_CENUMBER_HASH "$CEnumber\n";
    }
    elsif ($_ =~ /^>.+locus:/) 
    {
        print "Failed to scan locus from input line:\n";
        print "$_";
        print "Quitting.\n";
        die;
    }
    elsif ($_ =~ /^>(\S+)\s+(CE\d+)\D{1}/)
    {
        $cds = $1;
        $CEnumber = $2;
        $cds_hash_friendly = $cds;
        $cds_hash_friendly =~ s/\./\_/g;
        print CDS_TO_LOCUS_HASH    "$cds_hash_friendly\n";
        print CDS_TO_LOCUS_HASH    "$cds\n";
        print CDS_TO_CENUMBER_HASH "$cds_hash_friendly\n";
        print CDS_TO_CENUMBER_HASH "$CEnumber\n";
    }
    elsif ($_ =~ /^>/)
    {
        print "WARNING: The input line \"$_\" could not be interpreted by this script.\n";
    }
}
