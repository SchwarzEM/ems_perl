#!/usr/bin/perl -w

# clean_up_PfalPep.pl
# Erich Schwarz <emsch@its.caltech.edu>, 11/27/02.

# Purpose: Reformat PfalPep headers so that simple protein names are clearly used.
# PfalPep == Pfa3D7_WholeGenome_Annotated_PEP_2002.10.03-v2.fasta 09-Oct-2002 08:05   4.6M  [not 15.1M version!]
# from       http://www.plasmodb.org/restricted/data/P_falciparum/WG/cds.aa

# Sample header input:
#
# >Pfa3D7|pfal_chr1|PFA0005w|Annotation|Sanger (protein coding) erythrocyte membrane protein 1 (PfEMP1) Location=join(29733..34985,36111..37349)
#
# Sample desired header output:
# 
# >PFA0005w    Pfa3D7|pfal_chr1|PFA0005w|Annotation|Sanger (protein coding) erythrocyte membrane protein 1 (PfEMP1) Location=join(29733..34985,36111..37349)

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = $input . ".reformatted_headers";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    chomp ($input_line = $_);
    if ($input_line =~ />([^|]*)\|([^|]*)\|([^|]*)\|(.*)/)
    {
        print OUTPUT ">$3    $1\|$2\|$3\|$4\n";
    }
    else 
    {
        print OUTPUT "$input_line\n";
    }
}

close INPUT;
close OUTPUT;
