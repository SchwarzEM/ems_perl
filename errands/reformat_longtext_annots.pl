#!/usr/bin/perl -w

use strict;

# reformat_longtext_annots.pl
# by Erich Schwarz <emsch@its.caltech.edu>, 6/17/04

# Purpose: reformat human-writeable but obsolete annotations into newer .ace format (machine-better but human-unreadable).

# Latest bugfix: reformat "Gene : " annotations.

my $input                 = "";
my $output                = "";
my $reading_in_references = "no";
my $reading_in_longtext   = "no";
my $input_line            = "";
my $longtext_line         = "";
my @longtext_array        = "";
my $i                     = "";
my $annotation_textline   = "";
my $descript_type         = "";
my $old_longtext_line     = "";
my @paper_references      = "";
my $reference_to_print    = "";
my $indiv_reference       = "";


if ($#ARGV != 0) 
{
    print "LongText .ace file to be reformatted?: ";
    chomp($input = <STDIN>);
    $output = $input . ".reformatted_dot_ace";
} 
else 
{
    chomp($input = $ARGV[0]);
    $output = $input . ".reformatted_dot_ace";
}

open (INPUT, "$input") || die "Couldn't open input file. $!\n";
open (OUTPUT, ">$output") || die "Couldn't open output file. $!\n";

while (<INPUT>) 
{
    chomp ($input_line = $_);

    if ($input_line =~ /LongTextEnd/)
    {
        $reading_in_longtext = "no";

        # Time to rework a long, long annotation text line (>=1 sentences).
        # Clean up extra spaces; protect species names and oddities like "1st ed." and "et al.".

        $longtext_line =~ s/\s+/ /g;
        $longtext_line =~ s/\"/\'/g;    # added to keep double-quotes out of annotation text
        $longtext_line =~ s/^\s+//;
        $longtext_line =~ s/C\. elegans/C__elegans/g;
        $longtext_line =~ s/C\. briggsae/C__briggsae/g;
        $longtext_line =~ s/D\. discoideum/D__discoideum/g;
        $longtext_line =~ s/H\. sapiens/H__sapiens/g;
        $longtext_line =~ s/D\. melanogaster/D__melanogaster/g;
        $longtext_line =~ s/S\. cerevisiae/S__cerevisiae/g;
        $longtext_line =~ s/S\. pombe/S__pombe/g;
        $longtext_line =~ s/A\. nidulans/A__nidulans/g;
        $longtext_line =~ s/E\. coli/E__coli/g;
        $longtext_line =~ s/B\. subtilis/B__subtilis/g;
        $longtext_line =~ s/A\. suum/A__suum/g;
        $longtext_line =~ s/A\. thaliana/A__thaliana/g;
        $longtext_line =~ s/1st\. ed/1st__ed/g;
        $longtext_line =~ s/d\. ed/d__ed/g;
        $longtext_line =~ s/et al\./et al__/g;
        $longtext_line =~ s/e\. g\./e__g__/g;
        $longtext_line =~ s/Fig\. /Fig__/g;
        $longtext_line =~ s/deg\. C /degrees C /g;
        $longtext_line =~ s/deg\. C\./degrees C\./g;
        $longtext_line =~ s/[\.]{2,}/\./g;

        # Having de-periodized species names, split into an annotation sentence array.

        @longtext_array = split /\. /, $longtext_line;
     
        $i = 0;   # (re-)initialize "$i" for annotation sentence array.

        if ($longtext_array[0] eq "") 
        {
            print "Failed to detect longtext as sentence array.\n";

            # This is a pretty important check against formatting errors.
        }
        else 
        {
            until ($i > $#longtext_array)
            { 

                # Reformat back from protected periods (e.g., species names) to normal writing.

                $longtext_array[$i] =~ s/C__elegans/C\. elegans/g;
                $longtext_array[$i] =~ s/C__briggsae/C\. briggsae/g;
                $longtext_array[$i] =~ s/D__discoideum/D\. discoideum/g;
                $longtext_array[$i] =~ s/H__sapiens/H\. sapiens/g;
                $longtext_array[$i] =~ s/D__melanogaster/D\. melanogaster/g;
                $longtext_array[$i] =~ s/S__cerevisiae/S\. cerevisiae/g;
                $longtext_array[$i] =~ s/S__pombe/S\. pombe/g;
                $longtext_array[$i] =~ s/A__nidulans/A\. nidulans/g;
                $longtext_array[$i] =~ s/E__coli/E\. coli/g;
                $longtext_array[$i] =~ s/B__subtilis/B\. subtilis/g;
                $longtext_array[$i] =~ s/A__suum/A\. suum/g;
                $longtext_array[$i] =~ s/A__thaliana/A\. thaliana/g;
                $longtext_array[$i] =~ s/1st__ed/1st\. ed/g;
                $longtext_array[$i] =~ s/d__ed/d\. ed/g;
                $longtext_array[$i] =~ s/et al__/et al\. /g;
                $longtext_array[$i] =~ s/e__g__/e\. g\./g;
                $longtext_array[$i] =~ s/Fig__/Fig\. /g;
                $longtext_array[$i] =~ s/[\.]$//;

                $annotation_textline = $descript_type . " " . "\"$longtext_array[$i]\.\"";

                if ($i == 0)
                {
                    print OUTPUT "Concise_description";
                    print OUTPUT " ";
                    print OUTPUT "\"$longtext_array[$i]\.\"\n";
                }

                elsif ($i > 0) 
                {
                    print OUTPUT "$annotation_textline\n";
                }

                $i += 1;

                if ($i == 1)    # This means that the first line alone gets all the references.
                                # This was a quick-and-dirty solution.  It needs to be replaced by a better 
                                #   solution involving individually referenced annotation sentences.
                {
                    foreach $reference_to_print (@paper_references)
                    {
                        unless ($reference_to_print eq "")
                        {
                            print OUTPUT "$descript_type" . " " . "\"$longtext_array[0]\.\" ";
                            print OUTPUT "$reference_to_print\n";
                        }
                    }
                }
            }
        }
    }

    elsif ($reading_in_longtext eq "yes")
    {
        $old_longtext_line = $longtext_line;
        $longtext_line = $old_longtext_line . " " . $input_line;
    }

    elsif ($input_line =~ /^Gene :/ 
            || $input_line =~ /^Locus :/
            || $input_line =~ /^CDS :/ 
            || $input_line =~ /^Transcript :/ 
            || $input_line =~ /^Sequence :/) 
    {
        print OUTPUT "\n";
        print OUTPUT "$input_line\n";

        $descript_type         = "";
        @paper_references      = "";
        $reading_in_references = "yes";
        $reading_in_longtext   = "no";
        $longtext_line         = "";
        @longtext_array        = "";

    }

# needs to cope with stuff like:
# Detailed_description "chi-2" PMID_evidence "11163442"

    elsif ( ($input_line =~ /^Provisional_description\s+\"\S+\"\s+(\S+_evidence.+)/) && ($reading_in_references eq "yes") )
    {
        $indiv_reference = $1;   
        chomp ($indiv_reference);
        $descript_type = "Provisional_description";
        push (@paper_references, $indiv_reference);
    }

    elsif ( ($input_line =~ /^Detailed_description\s+\"\S+\"\s+(\S+_evidence.+)/) && ($reading_in_references eq "yes") )
    {
        $indiv_reference = $1;
        chomp ($indiv_reference);
        $descript_type = "Detailed_description";
        push (@paper_references, $indiv_reference);
    }

    elsif ($input_line =~ /^LongText :/) 
    {
        $reading_in_references = "no";
        $reading_in_longtext = "yes";
    }
    elsif ($input_line =~ /^LongText:/)  # quite common typographical error, so this error-check is important
    {
        print "Warning!  There is an improperly formatted \'LongText:\' entry in the .ace record!\n";
        die;
    }
}

close INPUT;
close OUTPUT;
