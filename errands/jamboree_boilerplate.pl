#!/usr/bin/perl -w

# jamboree_boilerplate.pl
# Purpose: generate very many frames in which to do jamboree annotations.
# Erich Schwarz <emsch@its.caltech.edu>, 1/29/03.

# Sample input:

# xyz-1

# Sample output:

# [empty line]
# Locus : "xyx-1"
# Provisional_description "xyz-1" Person_evidence "Lucky-Curator A"   // MANDATORY so that I can keep track of this!
# // Provisional_description "xyz-1" PMID_evidence "11111111"    // if you have a PMID reference you want to cite
# // Provisional_description "xyz-1" Paper_evidence "[cgc1111]"  // if you have a cgc paper you want to cite
# 
# LongText : "xyz-1"
# 
# Mindless jamboree boilerplate annotation sentence.
# Second mindless jamboree boilerplate annotation sentence.
# 
# ***LongTextEnd***
# [empty line]

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = $input . ".jamboree-proto.ace";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    print OUTPUT "\n";
    print OUTPUT "Locus : \"$input_line\"\n";
    print OUTPUT "Provisional_description \"$input_line\" Person_evidence \"Schwarz EM\"\n";
    print OUTPUT "\/\/ Provisional_description \"$input_line\" PMID_evidence \"11111111\"\n";
    print OUTPUT "\/\/ Provisional_description \"$input_line\" Paper_evidence \"[cgc1111]\"\n";
    print OUTPUT "\n";
    print OUTPUT "LongText : \"$input_line\"\n";
    print OUTPUT "\n";
    print OUTPUT "Mindless jamboree boilerplate annotation sentence.\n";
    print OUTPUT "Second mindless jamboree boilerplate annotation sentence. \n";
    print OUTPUT "\n";
    print OUTPUT "\*\*\*LongTextEnd\*\*\*\n";
    print OUTPUT "\n";
}

close INPUT;
close OUTPUT;
