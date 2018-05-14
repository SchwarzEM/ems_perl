#!/usr/bin/perl -w

# extract_m-6_blast_hits.pl
# Erich Schwarz, 5/9/02.
# Designed to extract non-redundant identifiers out of Blast searches
#   that did indeed give some sort of actual hit.

# Datestamping is used to ensure that a given analyses' files
#   are kept sorted out; two different runs of the script will
#   automatically generate distinctly stamped files and directories.

print "\n";
print "    Please enter a string that reliably defines\n";
print "    a set of Blast outputs in this directory, in \"-m 6\" format,\n";
print "    to be scanned for non-redundant hit identifiers.\n";
print "    E.g. \"ensembl.genscan.fa_5e-02.align-6\".\n";
print "\n";
print "    The string is?:   ";

$defining_string=<STDIN>;
chomp($defining_string);

print "\n";
print "    Please enter a text string briefly *labelling* the\n";
print "    Blast outputs being scanned. It should not match the\n";
print "    definition string.\n";
print "\n";
print "    The text string is?:   ";

$label_string=<STDIN>;
chomp($label_string);

if ($label_string =~ /$defining_string/) 
{
    die "Label string matched defining string!  Try again!\n";
}

$date = `date +"%s"`;
chomp($date);

$filename = "";

$rough_gi_list = $date . "." . $label_string . ".rough_gi_list";
$rough_id_list = $date . "." . $label_string . ".rough_id_list";

$gi_list = $date . "." . $label_string . ".gi_list";
$id_list = $date . "." . $label_string . ".id_list";

open (ROUGH_GI_LIST, ">$rough_gi_list");
open (ROUGH_ID_LIST, ">$rough_id_list");


# [Format of things to extract:]
#                                                                    Score     E
# Sequences producing significant alignments:                        (bits)  Value
# 
# JGI_5253                                                               42  7e-04
# 
# QUERY    55  LMRIVPPLAALILFCTYVLPLWGSGPQWNLVVGHHADICKKNWWRNLLFIHNYFGFSEM- 113
# JGI_5253 282 LLGIQPLHVFITLLTTGIVSLVQWGPYWFPFVDIIMD-CKSYWWANVVLVSNIIPVQQII 340
# [space, etc.]

# or:

#                                                                    Score     E
# Sequences producing significant alignments:                        (bits)  Value
# 
# gi|3805492|gb|AI223289.1|AI223289 qg53g06.x1 Soares_testis_NHT H...   173  3e-42
# gi|14455739|gb|BI049117.1|BI049117 PM2-UM0053-200401-018-b12 UM0...    42  0.007
# 
# QUERY    204 SSPGLVNRVLSWDIWSFLSSISYARYLVHPILIILYNGLQETLIHHTDTNMFYLFSGHRV 263
# 3805492  420 ShgGLVNRVLSWDIWSFLSSISYARYLVHPILIILYNGLQETLIHHTDTNMFYLFSGHRV 241
# 14455739 30  -----------------------------------------------------------V 32
# [space, etc.]

$scan_hits = "no";   # start reading hits after "Sequences producing..."
                       #  stop reading hits after QUERY

foreach (`ls | grep $defining_string`)
{
    $filename = $_;
    chomp($filename);
    open (FILENAME, "$filename");
    while (<FILENAME>)
    {
        if ($_ =~ /^Sequences producing significant alignments/)
        {
            $scan_hits = "yes";
        }
        elsif ($_ =~ /QUERY/) 
        {
            $scan_hits = "no";
        }
        elsif ($_ =~ /^gi\|(\d+)\|/) 
        {
            $gi_number = $1;
            if ($scan_hits eq "yes") 
            {
                print ROUGH_GI_LIST "$gi_number\n";
            }
        }
        elsif ($_ =~ /^(\S+)\s+/) 
        {
            $id_number = $1;
            if ($scan_hits eq "yes") 
            {
                print ROUGH_ID_LIST "$id_number\n";
            }
        }
    }
    print "The file $filename was scanned.\n";
}

close ROUGH_GI_LIST;
close ROUGH_ID_LIST;

# Tidy up the first output files.

system "sort $rough_gi_list | uniq - > $gi_list";
system "sort $rough_id_list | uniq - > $id_list";

system "rm $rough_gi_list";
system "rm $rough_id_list";

print "\n";
print "Job done.\n";
print "\n";
