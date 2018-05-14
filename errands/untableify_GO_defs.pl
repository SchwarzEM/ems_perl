#!/usr/bin/perl -w

# untableify_GO_defs.pl
# Erich Schwarz, 1/27/02.
# To un-table-ify extracted terms that I've worked on.

# GO.defs format:

# term: 1-phosphatidylinositol-4-phosphate kinase, class IB
# goid: GO:0004433
# definition: A class I PI3K activated via heterotrimeric G-proteins.
# definition_reference: MEDLINE:20504572


print "GO.defs file to untable-ify? ";
$input = <STDIN>;
chomp ($input);
$output = $input . ".untable_ified";

open (INPUT, "$input") || die "Can't open $input. $!";
open (OUTPUT, ">$output") || die "Can't open $output. $!";

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /\t/) 
    {
        @lines_to_print = split /\t/, $input_line;
        $i = 0;
        while ($i < $#lines_to_print) 
        {
            print OUTPUT "$lines_to_print[$i]";
            print OUTPUT "\n";
            $i++;
        }
    }
    else 
    {
        print OUTPUT "$input_line";
        print OUTPUT "\n";
    }
}
