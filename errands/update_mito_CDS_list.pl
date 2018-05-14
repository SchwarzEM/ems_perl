#!/usr/bin/perl -w

# update_mito_CDS_list.pl
# Erich Schwarz, 1/31/02
# Purpose: update old CDS names to wormpep72+ names.

print "Input file?: ";
$input_CDS_list = <STDIN>;
chomp ($input_CDS_list);

print "Input wormpepN?: ";
$input_wormpep = <STDIN>;
chomp ($input_wormpep);

$output_list = $input_CDS_list . ".output";
$rough_output_list = $input_CDS_list . ".rough_output";

open (INPUT_CDS, "$input_CDS_list");
open (ROUGH_OUTPUT_CDS, ">$rough_output_list");

while (<INPUT_CDS>) 
{
    foreach ($_) 
    {
        $input_cds = $_;
        chomp ($input_cds);
        open (INPUT_WPEP, "$input_wormpep");
        while (<INPUT_WPEP>) 
        {
            foreach ($_) 
            {
                $input_wpep_line = $_;
                if ($input_wpep_line =~ /^>(\S+\d+)\s+/) 
                {
                    $input_wpep_cds = $1;
                    if ($input_cds eq $input_wpep_cds) 
                    {
                        print ROUGH_OUTPUT_CDS "$input_cds\n";
                    }
                }
                 elsif ($input_wpep_line =~ /^>(\S+\d+[A-Z]{1})\s+/) 
                {
                    $input_wpep_cds = $1;
                    $input_wpep_cds_trimmed = $input_wpep_cds;
                    $input_wpep_cds_trimmed =~ s/^(\S+\d+)[A-Z]{1}/$1/;
                    if ($input_cds eq $input_wpep_cds_trimmed) 
                    {
                        print ROUGH_OUTPUT_CDS "$input_wpep_cds\n";
                    }
                }
            }
        }
    }
}

close INPUT_CDS;
close INPUT_WPEP;
close ROUGH_OUTPUT_CDS;

system "sort $rough_output_list | uniq - > $output_list";
system "rm $rough_output_list";
