#!/usr/bin/perl

# Program: wormpep2go_to_locus_list.pl
#    Erich Schwarz, 6/29/01
# 
# Purpose: Stupid fast hack to get locus list.

# Load files.

# Default hashes to load:
#     component.ontology.hash  
#     function.ontology.hash
#     process.ontology.hash
#     wormpep.CDS-to-locus.hash
# These files are all one-line-per-item lists.

print "Please give your name and e-mail address: ";
$authoraddress = <STDIN>;
chomp($authoraddress);

print "What will input wormpep2go file be?\n";
print "Default is \"wormpep2go\": ";
$infile = <STDIN>;
chomp ($infile);
unless ($infile =~ /\w+/)
{
    $infile = "wormpep2go";
}
$outfile = ("locuslist_for_gene_association.wb");

open (INFILE, "$infile") || die "wormpep2go file $infile not found\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outfile\n";

print "\n";
print "The input file is $infile; the output file $outfile\n";
print "\n";

print "What will input component.ontology.hash file be?\n";
print "Default is \"component.ontology.hash\": ";
$comp_file = <STDIN>;
chomp ($comp_file);
unless ($comp_file =~ /\w+/)
{
    $comp_file = "component.ontology.hash";
}

print "What will input function.ontology.hash file be?\n";
print "Default is \"function.ontology.hash\": ";
$funct_file = <STDIN>;
chomp ($funct_file);
unless ($funct_file =~ /\w+/)
{
    $funct_file = "function.ontology.hash";
}

print "What will input process.ontology.hash file be?\n";
print "Default is \"process.ontology.hash\": ";
$proc_file = <STDIN>;
chomp ($proc_file);
unless ($proc_file =~ /\w+/)
{
    $proc_file = "process.ontology.hash";
}

print "What will input wormpep.CDS-to-locus.hash be?\n";
print "Default is \"wormpep.CDS-to-locus.hash\": ";
$cds_to_locus_file = <STDIN>;
chomp ($cds_to_locus_file);
unless ($cds_to_locus_file =~ /\w+/)
{
    $cds_to_locus_file = "wormpep.CDS-to-locus.hash";
}

system "touch ontology.pre-hash";
system "cat $comp_file >> ontology.pre-hash";
system "cat $funct_file >> ontology.pre-hash";
system "cat $proc_file >> ontology.pre-hash";

open (ONTOLOGY_PRE_LIST, "ontology.pre-hash") || die "can't open ontology.pre-hash\n";
open (ONTOLOGY_LIST, ">ontology.hash") || die "can't open ontology.hash\n";

while (<ONTOLOGY_PRE_LIST>) 
    {
        $ontology_term = $_;
        chomp($ontology_term);
        $ontology_term =~ s/\:/\_/g;
        print ONTOLOGY_LIST "$ontology_term\n";
    }

close ONTOLOGY_PRE_LIST;
close ONTOLOGY_LIST;

system "rm ontology.pre-hash";

open (GO_TERM_TO_ONTOLOGY_LIST, "ontology.hash") || die "can't open ontology.hash\n";
open (CDS_TO_LOCUS_LIST, "$cds_to_locus_file") ||  die "can't open $cds_to_locus_file\n";

print "\n";
print "The ontology hash file is \"ontology.hash\"\n";
print "The wormpep.CDS-to-locus.hash is \"$cds_to_locus_file\"\n";
print "\n";

# 1. In wormpep2go:
#
#     First column is CDS.  Save.
#         Use as gene name if locus not assigned.
#     Second column is WP number -- save, use as database ID number.
#     Third column is Interpro entry: keep (to modify evidence code).
#     Fourth column is GO number: keep.
#
# 2. Spit out, for each line of wormpep2go:
#
#     Column 1: WB  [Invariant database abbreviation.]
#     Column 2: WP:---  [ID number, from col. 2 of wormpep2go.]
#     Column 3: Output of CDS-to-gene hash, given CDS.  [Gene.]
#     Column 4: GO term, taken straight.
#     Column 5: PUBMED:11159333  [Most recent Interpro reference.]
#     Column 6: IEA  [Evidence -- Inferred by Electronic Annotation.]
#     Column 7: Ontology letter (P, F, C): from GO hash lookup.
#     Column 8: CDS name, and gene name if it's distinct, separated by stave.
#                  Eventually we want classical synonyms in here too.

@go_term_to_ontology_list = <GO_TERM_TO_ONTOLOGY_LIST>;
@cds_to_locus_list = <CDS_TO_LOCUS_LIST>;

chomp(@go_term_to_ontology_list);
chomp(@cds_to_locus_list);

%go_term_to_ontology_hash = @go_term_to_ontology_list;
%cds_to_locus_hash = @cds_to_locus_list;

# Processing of INFILE must cope with this sort of input:
# 
# 2L52.1	WP:CE20433	INTERPRO:IPR000822	GO:0003700
# 2L52.1	WP:CE20433	INTERPRO:IPR000822	GO:0005634
# 2L52.1	WP:CE20433	INTERPRO:IPR000822	GO:0006355
# 4R79.1	WP:CE19649	INTERPRO:IPR000130	GO:0006508

# Note: perl hashes appear to choke on strings containing "." or ":".

# Remaining bug:
# Use of uninitialized value in string at ./wormpep2go_to_GO-gene-assignments.pl line 178, <INFILE> line 20966.
# Use of uninitialized value in concatenation (.) at ./wormpep2go_to_GO-gene-assignments.pl line 189, <INFILE> line 20966.

while (<INFILE>) 
{
    if ($_ =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)/) 
    {
        $cds_literal = $1;
        $wp_id = $2;
        $interpro_code = $3;
        $go_term_literal = $4;
        $cds_hash_friendly = $cds_literal;
        $cds_hash_friendly =~ s/\./\_/g;
        $go_term_hash_friendly = $go_term_literal;
        $go_term_hash_friendly =~ s/\:/\_/g;

        # Print "Gene_symbol" (1st) column.
        print OUTFILE "$cds_to_locus_hash{$cds_hash_friendly}";
        print OUTFILE "\n";
    }
}
