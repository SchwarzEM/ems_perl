#!/usr/bin/perl

# GOify_worm_prions.pl
# 
# Erich Schwarz, 1/30/02
# 
# Purpose: make GO-format files from a wormpepN DIANA hitlist.

# #################### run with standard defaults ####################

$defaults = "yes";

print "\n";
print "Run with standard defaults?  [Default is \"yes\"]: ";
$defaults = <STDIN>;
$defaults =~ tr/A-Z/a-z/;
unless ($defaults =~ /no/) 
{
    print "\n";
    print "OK, running with standard defaults...\n";

    $prion_hitlist        = "diana-hitlist";	# was: $seq_list
    $cds_to_locus_file    = "wormpep.CDS-to-locus.hash";
    $cds_to_CEnumber_file = "wormpep.CDS-to-CEnumber.hash";

    $gene_assoc_outfile   = "prion_gene_association.wb";
}

# #################### run with hand-entered defaults ####################

elsif ($defaults =~ /no/) 
{
    print "\n";
    print "Please name the DIANA prion-like hitlist to use (one line per sequence).\n";
    print "The default is \"diana-hitlist\": ";
    $prion_hitlist=<STDIN>;
    chomp($prion_hitlist);
    unless ($prion_hitlist =~ /\w+/)
    {
        $prion_hitlist = "diana-hitlist";
    }

    print "\n";
    print "Please name the wormpep.CDS-to-locus.hash to use (default is \"wormpep.CDS-to-locus.hash\"): ";
    $cds_to_locus_file = <STDIN>;
    chomp ($cds_to_locus_file);
    unless ($cds_to_locus_file =~ /\w+/) 
    {
        $cds_to_locus_file = "wormpep.CDS-to-locus.hash";
    }

    print "\n";
    print "Please name the wormpep.CDS-to-CEnumber.hash to use (default is \"wormpep.CDS-to-CEnumber.hash\"): ";
    $cds_to_CEnumber_file = <STDIN>;
    chomp ($cds_to_CEnumber_file);
    unless ($cds_to_CEnumber_file =~ /\w+/)
    {
        $cds_to_CEnumber_file = "wormpep.CDS-to-CEnumber.hash";
    }

    print "\n";
    print "Please give the name that the prion gene association output file should have;\n";
    print "The default choice would be \"prion_gene_association.wb\": ";
    $gene_assoc_outfile = <STDIN>;
    chomp($gene_assoc_outfile);
    unless ($gene_assoc_outfile =~ /\w+/)
    {
        $gene_assoc_outfile = "prion_gene_association.wb";
    }
}

# #################### abort run if defaults not entered ####################

else 
{
    print "Hm.  I didn't get \"yes\" or \"no\" for the defaults.  Try again?\n";
    die;
}

# #################### announce working defaults ####################

print "\n";
print "The DIANA prion-like hitlist will be:        $prion_hitlist\n";
print "The wormpep.CDS-to-locus.hash is:            $cds_to_locus_file\n";
print "The wormpep.CDS-to-CEnumber.hash is:         $cds_to_CEnumber_file\n";
print "\n";
print "The gene association output file will be:    $gene_assoc_outfile\n";
print "\n";

# set up sortable rough output file

$rough_gene_assoc_outfile = "rough_" . $gene_assoc_outfile;

# #################### make the date-stamp for GO files ####################

$datestamp = `date -I`;
chomp ($datestamp);
$orig_datestamp = $datestamp;
$datestamp =~ s/-//g;

# #################### get Interpro-derived gene associations for GO ####################

open (PRION_HITLIST, "$prion_hitlist") 
    || die "Wormpep2go file $worm2go_infile not found. $!\n";
open (CDS_TO_LOCUS_LIST, "$cds_to_locus_file") 
    || die "Can't open $cds_to_locus_file -- $!\n";
open (CDS_TO_CENUMBER_LIST, "$cds_to_CEnumber_file") 
    || die "Can't open $cds_to_CEnumber_file -- $!\n";

open (ROUGH_GENE_ASSOC, ">$rough_gene_assoc_outfile") 
    || die "Couldn't open file $rough_gene_assoc_outfile. $!\n";

# 2. Spit out:
#
#     Column  1: WB  [Invariant database abbreviation.]
#     Column  2: WP:---  [ID number, from col. 2 of wormpep2go.]
#     Column  3: Output of CDS-to-gene hash, given CDS.  [Gene.]
#     Column  4: [NOT] column!    <-- must be sure this is in.
#     Column  5: GO term: for this it's automatically "GO:xxx == epig"; the program reminds this.
#     Column  6: PUBMED:11050225  [I.e., Mich. and Weissman DIANA reference.]
#     Column  7: ISS  [Evidence -- Inferred by Similarity or Structure.]
#     Column  8: [With] -- used for various things like IGI, IPI. <-- mandatory skip.
#     Column  9: Ontology letter (P, F, C): for this it's automatically "P".
#     Column 10: [Name(|Name)]: in most cases this is pro forma, thus blank.
#     Column 11: [Synonym(|Synonym)], if one or more such exists.
#     Column 12: [DB_Object_Type]: in this case necessarily "protein".
#     Column 13: [Taxon_ID]: always "taxon:6239".
#

@cds_to_locus_list = <CDS_TO_LOCUS_LIST>;
@cds_to_CEnumber_list = <CDS_TO_CENUMBER_LIST>;

chomp(@cds_to_locus_list);
chomp(@cds_to_CEnumber_list);

%cds_to_locus_hash    = @cds_to_locus_list;
%cds_to_CEnumber_hash = @cds_to_CEnumber_list;

# Note: perl hashes appear to choke on strings containing "." or ":".

while (<PRION_HITLIST>) 
{
    $cds_literal = $_;
    chomp ($cds_literal);
    $cds_hash_friendly = $cds_literal;
    $cds_hash_friendly =~ s/\./\_/g;

    # Print "DB" (1st) column.
    print ROUGH_GENE_ASSOC "WB";
    print ROUGH_GENE_ASSOC "\t";

    # Print "DB_Object_ID" (2cd) column -- i.e., "Gene_ID".
    print ROUGH_GENE_ASSOC "$cds_to_CEnumber_hash{$cds_hash_friendly}";
    print ROUGH_GENE_ASSOC "\t";

    # Print "DB_Object_Symbol" (3rd) column -- i.e., "Gene_Symbol".
    print ROUGH_GENE_ASSOC "$cds_to_locus_hash{$cds_hash_friendly}";
    print ROUGH_GENE_ASSOC "\t";

    # Print "NOT" (4th) column.
    # In most cases this is pro forma, thus blank.
    print ROUGH_GENE_ASSOC "\t";

    # Print "GOid" (5th) column.
    print ROUGH_GENE_ASSOC "GO:0040030";         # Hard-wired for this particular output:
    print ROUGH_GENE_ASSOC "\t";                 # "epigenetic control of protein activity ; GO:0040030"

    # Print "DB:Ref" (6th) column.
    print ROUGH_GENE_ASSOC "PUBMED:11050225";    # Hard-wired for this particular output: 
    print ROUGH_GENE_ASSOC "\t";                 # PUBMED:11050225  -- Michelitsch and Weissman (2000), PNAS 97, 11910-11915.

    # Print "Evidence" (7th) column.
    print ROUGH_GENE_ASSOC "ISS";
    print ROUGH_GENE_ASSOC "\t";

    # Print "With" (8th) column.
    # In most cases this is pro forma, thus blank.
    print ROUGH_GENE_ASSOC "\t";

    # Print "Aspect" (9th) column.
    print ROUGH_GENE_ASSOC "P";    # In this case, it's an automatic "P".
    print ROUGH_GENE_ASSOC "\t";

    # Print "Name(|Name)" (10th) column.
    # In most cases this is pro forma, thus blank.
    print ROUGH_GENE_ASSOC "\t";

    # Print "Synonym(|Synonym)" (11th) column.
    if ($cds_literal eq $cds_to_locus_hash{$cds_hash_friendly}) 
    {
        print ROUGH_GENE_ASSOC "\t";
    }
    else 
    {
        print ROUGH_GENE_ASSOC "$cds_literal";
        print ROUGH_GENE_ASSOC "\t";
    }
    # Print "DB_Object_Type" (12th) column -- usually "protein".
    # But note that, for RNA products of genes, one needs "transcript".
    print ROUGH_GENE_ASSOC "protein";
    print ROUGH_GENE_ASSOC "\t";

    # Print "Taxon_ID" (13th) column -- always "taxon:6239".
    print ROUGH_GENE_ASSOC "taxon:6239";
    print ROUGH_GENE_ASSOC "\t";

    # Print "Date" (14th) column -- e.g., "20020127" for Jan. 27, 2002.
    print ROUGH_GENE_ASSOC "$datestamp";
    print ROUGH_GENE_ASSOC "\n";
}

close PRION_HITLIST;
close CDS_TO_LOCUS_LIST;
close CDS_TO_CENUMBER_LIST;
close GENE_ASSOC;

system "sort $rough_gene_assoc_outfile | uniq - > $gene_assoc_outfile";
system "rm $rough_gene_assoc_outfile";
