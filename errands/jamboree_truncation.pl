#!/usr/bin/perl -w

# jamboree_truncation.pl
# Purpose: give me an easy way to chop off already-done and half-done jamboree annotations.
# Erich Schwarz <emsch@its.caltech.edu>, 2/11/03.

print "input file: ";
$input = <STDIN>;
chomp($input);
$output = $input . ".truncated_jambo-annots";

$should_I_keep_printing = "no";

open (INPUT, "$input") || die;
open (OUTPUT, ">$output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /stop printing here/i) 
    {
        $should_I_keep_printing = "no";
    }
    elsif ($input_line =~ /start printing here/i) 
    {
        $should_I_keep_printing = "yes";
    }
    elsif ($should_I_keep_printing eq "yes") 
    {
        print OUTPUT "$input_line\n";
    }
}

close INPUT;
close OUTPUT;
