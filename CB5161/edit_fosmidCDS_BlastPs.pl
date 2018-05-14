#!/usr/bin/perl

# edit_fosmidCDS_BlastPs.pl
# Erich Schwarz <emsch@its.caltech.edu>, 10/26/03
# Purpose: get reasonable format of best hits from BlastP/wormpep scans of predicted fosmid proteins

if ($#ARGV != 0) 
{
    print "BlastP output to edit? ";
    chomp ($input_file = <STDIN>);
} 
else 
{
    $input_file = $ARGV[0];
}

$output_file = $input_file . ".extract";

open (INPUT, "$input_file") || die;
open (OUTPUT, ">$output_file") || die;

$seq_name = "";
$scan = "no";
$percent_id = "ND";

while (<INPUT>) 
{
    chomp($input_line = $_);

    if ($input_line =~ /^Query=\s+(\S+)\s+/)
    {
        unless ($aligned_res == 0) 
        {
            print OUTPUT "\tIdentities = $ident_res/$aligned_res ($percent_id%)\n";
        }
        if ($aligned_res == 0 && $seq_name != "")
        {
            print OUTPUT "\tIdentities = n/a";
        }
        unless ($seq_name == "") 
        {
            print OUTPUT "\n";
        }
        $scan = "start";
        $ident_res   = 0;
        $aligned_res = 0; 
        $seq_name = $1;
        print OUTPUT "$seq_name\t";
    }
    elsif (($scan eq "start") && ($input_line =~ /\s+\((\d+) letters\)/))
    {
        print OUTPUT "Length = $1\t";
    }

    elsif (($scan eq "start") && ($input_line =~ /^>(.+)\s+CE\d+\s+(locus:\S+)\s/))
    {
        print OUTPUT "$1 $2";
        $scan = "ongoing";
    }
    elsif (($scan eq "start") && ($input_line =~ /^>(.+)\s+CE\d+\s+/))
    {
        print OUTPUT "$1";
        print OUTPUT "             ";   # filler to square off lines
        $scan = "ongoing";
    }
    elsif (($scan eq "start") && ($input_line =~ /^>(.+) W =/))
    {
        print OUTPUT "$1";
        $scan = "ongoing";
    }
    elsif (($scan eq "start") && ($input_line =~ /^>(.+)$/)) 
    {
        print OUTPUT "$1";
        $scan = "ongoing";
    }

    elsif (($scan eq "ongoing") && ($input_line =~ /(Length = \d+)/))
    {
        print OUTPUT "\t$1";
    }
    elsif (($scan eq "ongoing") && ($input_line =~ /Identities = (\d+)\/(\d+)/))
    {
        $ident_res   = $ident_res   + $1;
        $aligned_res = $aligned_res + $2;  
        $percent_id  = int ( ($ident_res / $aligned_res) * 100)
    }

    elsif (($scan eq "ongoing") && ($input_line =~ /^>(.+)$/))
    {
        $scan = "no";
    }
}

unless ($aligned_res == 0)
{
    print OUTPUT "\tIdentities = $ident_res/$aligned_res ($percent_id%)";
}
    if ($aligned_res == 0 && $seq_name != "")
{
    print OUTPUT "\tIdentities = n/a";
}
print OUTPUT "\n";
