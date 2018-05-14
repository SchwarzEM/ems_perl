#!/usr/bin/perl -w

# Program: cgc_to_endnote.pl
#    Erich Schwarz, emsch@its.caltech.edu, 12/8/00
# 
# Purpose: Convert CGC references to Endnote references.

# 1. If not given a CGC file as argument, ask for its name.

if ($#ARGV != 0) 
{
	print "Required: input CGC file.\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} 
else 
{
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".endnote");
$botchfile = ($infile . ".botchlist");
print "The input file is $infile; the output file $outfile;\n";
print "   and the list of botched citations is $botchfile.\n";

# 2. Open the CGC and future Endnote files or die.

open (INFILE, "$infile") || die "CGC file $infile not found.  $!\n";
open (OUTFILE, ">$outfile") || die "Couldn't open Endnote file $outfile.  $!\n";
open (BOTCHFILE, ">$botchfile") || die "Couldn't open botch list $botchfile.  $!\n";

# 3. Insert header required for .endnote to be readable by Endnote
#       as a tab-delimited file.

print OUTFILE ("*Journal Article\n");
print OUTFILE ("Label\tAccession Number\tAuthor\tTitle\tJournal\tVolume\tPages\tYear\tAbstract\n");

# 4. Do a big initial round of initialization (to be redone at each Key).

$key = "";
$medline = "";  # Could alternatively be "no_medline_number".
$authors = "";
$title = "";
$authors_reading = "yes";	# Purpose: 
$title_reading = "no";          # Make sure that blocked 
$citation_reading = "no";       #   text is only read when 
$abstract_reading = "no";	#   it should be.

$volume_title = "";
$volume_number = "";
$page_numbers = "";
$year = "";
$abstract = "";
$first_citation_line = "";
$full_citation = "";
$append_to_citation = "";
$citation_text = "";

$botch_warning = "";

# 5. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
	if ($_ =~ /Key:[\D]+([\d]+)/) 
	{
		# Initialize a bunch of non-key variables.

		$medline = "";  # Could alternatively be "no_medline_number".
		$authors = "";
		$title = "";
		$authors_reading = "no";	# Purpose:
		$title_reading = "no";		# Make sure that blocked 
		$citation_reading = "no";	#   text is only read when 
                $abstract_reading = "no";	#   it should be.

		$volume_title = "";
		$volume_number = "";
		$page_numbers = "";  
		$year = "";
		$abstract = "";
		$first_citation_line = "";
		$full_citation = "";
		$append_to_citation = "";
		$citation_text = "";

		$botch_warning = "";

		# Get on with ordinary stuff.

		# Scan number after "Key: " into $key.
		$key = $1;
		if ($key eq "") 
		{
			die ("No CGC index key no. found! $!\n");
		}
	}
	elsif ($_ =~ /Medline:[\D]+([\d]+)/) 
	{
		# Scan number after "Medline: " into $medline.
		$medline = $1;
        }
	elsif ($_ =~ /Authors:[ ]+(.+)/) 
	{
		# Pick up list into $authors.
		$authors = $1;
		chomp ($authors);
		# Convert ";" into "//"
		$authors =~ s/;/\/\//g;  # Must be global subst.
		$authors =~ s/([a-zA-Z]+) ([a-zA-Z]+)/$1, $2/g;	#Global comma-ing.
		$authors_reading = "yes";
	}
        elsif ($_ =~ /^[ ]{13}(.+)/ && $authors_reading eq "yes")
        {
                # Append following lines of blocked text to $authors.
                $append_to_authors = $1;   
                chomp ($append_to_authors);
                $authors = $authors . " " . $append_to_authors;
                # Convert ";" into "//"
                $authors =~ s/;/\/\//g;  # Must be global subst.
                $authors =~ s/([a-zA-Z]+) ([a-zA-Z]+)/$1, $2/g; #Global comma-ing.
        }
        elsif ($_ =~ /Title: (.+)/)
        {
		$authors_reading = "no";
                # Enter all text after "Title: " into $title.
                $title = $1;
                chomp ($title);
		$title_reading = "yes";
	}
	elsif ($_ =~ /^[ ]{13}(.+)/ && $title_reading eq "yes") 
	{
		# Append following lines of blocked text to $title.		
		$append_to_title = $1;
		chomp ($append_to_title);
		$title = $title . " " . $append_to_title;
	}

# ---------------- Citation processing.  ----------------

        elsif ($_ =~ /Citation: (.+)/) 
	{
		$title_reading = "no";
		$citation_reading = "yes";
		$first_citation_line = $1;   
		chomp ($first_citation_line);
		$full_citation = $first_citation_line;
	}
        elsif ($_ =~ /^[ ]{13}(.+)/ && $citation_reading eq "yes") 
	{
                # Append following lines of blocked text to $full_citation.

                $append_to_citation = $1;   
                chomp ($append_to_citation);
                $full_citation = $full_citation . " " . $append_to_citation;
	}
        elsif ($citation_reading eq "yes")
        {
                # Process out values from $full_citation.

		if ($full_citation =~ /^(.+) : ([\d]+-[\d]+)[\s]+([\d]+)/) 
		{
			$citation_text = $1;
                	$page_numbers = $2;
                	$year = $3;
                	$citation_reading = "no";			
		}                
		elsif ($full_citation =~ /^(.+) ([\d]+): ([\d]+-[\d]+)[\s]+([\d]+)/)
        	{
                	$citation_text = $1;   
                	$volume_number = $2;  
                	$page_numbers = $3;   
                	$year = $4;
			$citation_reading = "no";
        	}
		else 
		{
		$botch_warning = "Citation $key was botched.\n";
		}
        }

# ---------------- End citation processing.  ----------------

# Skip "Type: " and "Genes: " lines and their text.

	elsif ($_ =~ /Abstract: (.+)/) 
	{
		# Enter all text after "Abstract: " into $abstract.

		$citation_reading = "no";
		$abstract = $1;
		chomp ($abstract);
		$abstract_reading = "yes";
	}
        elsif ($_ =~ /^[ ]{13}(.+)/ && $abstract_reading eq "yes")
        {
                # Append following lines of blocked text to $abstract.
                $append_to_abstract = $1;  
                chomp ($append_to_abstract);
                $abstract = $abstract . " " . $append_to_abstract;
        }

# 5. Print output.

	elsif ($abstract_reading eq "yes")
	{
		print OUTFILE ($key . "\t");
		print OUTFILE ($medline . "\t");
		print OUTFILE ($authors . "\t");
                print OUTFILE ($title . "\t");
		print OUTFILE ($citation_text . "\t");
		print OUTFILE ($volume_number . "\t");
		print OUTFILE ($page_numbers . "\t");
		print OUTFILE ($year . "\t");
		print OUTFILE ($abstract . "\n");
		if ($botch_warning ne "") 
		{
			print BOTCHFILE ($botch_warning . "\n");
		}
		$abstract_reading = "no";
		$botch_warning = "";
	}
}
