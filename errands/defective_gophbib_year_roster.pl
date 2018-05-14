#!/usr/bin/perl

use strict;

# gophbib_year_roster.pl
# Erich Schwarz <emsch@its.caltech.edu>, 10/17/02.
# Purpose: extract numbers of papers per year from gophbib.

# Example input line:  "   Citation: Nematologica 21: 151-162 1975"
# Desired output: "1975\tX\n", where X is total number of '1975' years in gophbib.

print "Input file: ";

chomp(my $input  = <STDIN>);
open (INPUT, "$input") || die;

my $rough_yearlist = $input . ".rough_yearlist";
open (ROUGH_YEARLIST, ">$rough_yearlist") || die;

while (<INPUT>) 
{
    chomp (my $input_line = $_);

    if ($input_line =~ /^   Citation: .+\s(\d+)$/)
    {
        print ROUGH_YEARLIST "$1\n";
    }   
}

close INPUT;
close ROUGH_YEARLIST;

my $yearlist = $input . ".yearlist";
system "sort $rough_yearlist > $yearlist";
open (YEARLIST, "$yearlist") || die;

my $last_year_read = "";
my $year_count = "";

my $output = $input . ".year_roster"; 
open (OUTPUT, ">$output") || die;

my $last_year_count = "";
my $absolute_increment = "";
my $percent_increment = "";

print OUTPUT "Year\t Papers\t +\/-\t +\/- %\n";

print OUTPUT "\n";

while (<YEARLIST>) 
{
    chomp (my $input_year = $_);
    if ($input_year != $last_year_read) 
    {
        unless ($last_year_read == "") 
        {
            $absolute_increment = $year_count - $last_year_count;
            if ($last_year_count != 0) 
            {
                $percent_increment = (($absolute_increment / $last_year_count) * 100);
            }
            printf OUTPUT "%4d\t%4d\t%4d\t%4.0f\n", 
                $last_year_read, $year_count, $absolute_increment, $percent_increment; 
            unless ($year_count == 0)
            {
                $last_year_count = $year_count;
            }
            if (($last_year_read % 10) == 0) 
            {
                print OUTPUT "\n";
            }
        }
        $last_year_read = $input_year;
        $year_count = 1;
    }
    elsif ($input_year == $last_year_read) 
    {
        $year_count += 1;
    }
}

printf OUTPUT "%4d\t%4d\n", $last_year_read, $year_count;

close YEARLIST;
close OUTPUT;

system "rm $rough_yearlist; rm $yearlist";
