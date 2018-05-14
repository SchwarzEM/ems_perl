#!/usr/bin/perl -w 

# clean_up_annot_PMIDrefs.pl
# Erich Schwarz <emsch@its.caltech.edu>, 4/28/03

# Purpose: convert annotation references from 'PMID_evidence "X"' to 'Paper_evidence "[pmidX]"'

print "Annotation file whose names are to be cleaned or warned about? ";

chomp ($input_file         = <STDIN>);
$rePMIDed_annot_file        = $input_file . ".rePMIDed_annot_file";

open (INPUT, "$input_file") || die;
open (REPMIDED_ANNOTS, ">$rePMIDed_annot_file") || die;

while (<INPUT>)
{
    chomp ($input_line = $_);
    if ($input_line =~ /(.+\s)PMID_evidence\s+\"(\d+)\"(.*)$/) 
    {
        $input_line_part_1 = $1;
        $input_line_part_2 = $2;
        $input_line_part_3 = $3;
        $input_line = $input_line_part_1 . "Paper_evidence \"\[pmid" . $input_line_part_2 . "\]\"" . $input_line_part_3;
    }
    print REPMIDED_ANNOTS "$input_line\n";
}

close INPUT;
close REPMIDED_ANNOTS;
