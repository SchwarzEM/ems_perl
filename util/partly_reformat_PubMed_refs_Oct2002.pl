#!/usr/bin/perl -w

# partly_reformat_PubMed_refs_Oct2002.pl  [permanent archival old version]
# Erich Schwarz <emsch@its.caltech.edu>, 10/20/02.
# Purpose: convert (e.g) "Lin JJ, Lobo NF, Lopez JR, Malek JA" to "Lin, J.J, Lobo, N.F., Lopez, J.R., Malek, J.A.".

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = $input . ".output";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    $input_text = $input_text . " " . $input_line;
}
close INPUT;

while ($input_text =~ /^([^,]+)\s([A-Z]+)\,(.+)/ )
{
    $surname   = $1;
    $initials  = $2;
    $left_text = $3;
    print OUTPUT "$surname, ";
    while ($initials =~ /([A-Z]{1})(.*)/) 
    {
            print OUTPUT "$1.";
            $initials = $2;
    }
    print OUTPUT ",";
    $input_text = $left_text;
}
print OUTPUT "$left_text\n";

close OUTPUT;
