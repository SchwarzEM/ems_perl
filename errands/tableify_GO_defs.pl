#!/usr/bin/perl -w

# tableify_GO_defs.pl
# Erich Schwarz, 1/27/02.
# To make it easy to extract terms I've worked on.

# GO.defs format:

# term: 1-phosphatidylinositol-4-phosphate kinase, class IB
# goid: GO:0004433
# definition: A class I PI3K activated via heterotrimeric G-proteins.
# definition_reference: MEDLINE:20504572


print "GO.defs file to table-ify? ";
$input = <STDIN>;
chomp ($input);
$output = $input . ".table_ified";

open (INPUT, "$input") || die "Can't open $input. $!";
open (OUTPUT, ">$output") || die "Can't open $output. $!";

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /^definition_reference/) 
    {
        print OUTPUT "$input_line";
        print OUTPUT "\n";
    }
    elsif ($input_line =~ /^\s+$/) 
    {
        print OUTPUT "\t";
    }
    else 
    {
        print OUTPUT "$input_line";
        print OUTPUT "\t";
    }
}

close INPUT;
close OUTPUT;
