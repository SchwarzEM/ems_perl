#!/usr/bin/perl -w

# check_funct_annot_members.pl
#
# Erich Schwarz <emsch@its.caltech.edu>, sometime after 4/17/02; latest revision 11/21/03.
# Purpose: get clean list of members of a functional annotation .ace file.

print "Input file: ";
$input = <STDIN>;
chomp($input);
$rough_output = $input . ".rough_list";
$output = $input . ".list";

open (INPUT, "$input") || die;
open (ROUGH_OUTPUT, ">$rough_output") || die;

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);

# example of two input lines:
# Locus : "mab-18"    // == 'F14F3.1a', 'F14F3.1b', 'F14F3.1c'
# Sequence : "B0261.2a"   // blah blah comment crud

    if ($input_line =~ /^(Locus : \"\S+\")/)
    {
        $output_line = $1;
        chomp ($output_line);
        print ROUGH_OUTPUT "$output_line\n";
    }
    elsif ($input_line =~ /^(Sequence : \"\S+\")/) 
    {
        $output_line = $1;
        chomp ($output_line);
        print ROUGH_OUTPUT "$output_line\n";
    }
    elsif ($input_line =~ /^(CDS : \"\S+\")/)
    {
        $output_line = $1;
        chomp ($output_line);
        print ROUGH_OUTPUT "$output_line\n";
    }
}

close INPUT;
close ROUGH_OUTPUT;

system "sort $rough_output > $output";
system "rm $rough_output";

print "The sorted list of entries in:  $input\n";
print "                           is:  $output\n";
