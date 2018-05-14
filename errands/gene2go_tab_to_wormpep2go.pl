#!/usr/bin/perl -w
#
# gene2go_tab_to_wormpep2go.pl
# Erich Schwarz, 10/11/01, with minor tweaks 1/25/02
#
# Two purposes:
# 
#   1. Process gene2go.tab from Tablemaker into 
#      digestable format for geneontology.org.
# 
#   2. Scan all wormpepN CDSes for INTERPRO motifs 
#      that *ought* to have hits but have not been 
#      correctly processed by the Sanger Centre.

print "What will input interpro2go file be?  [default: \"interpro2go\"] ";
$input_interpro2go = <STDIN>;
chomp ($input_interpro2go);
unless ($input_interpro2go =~ /\w+/)
{
    $input_interpro2go = "interpro2go";
}


open (INPUT_INTERPRO2GO, "$input_interpro2go") 
    || die "This program requires a readable interpro2go table! $!\n";
close INPUT_INTERPRO2GO;

print "What will input Tablemaker file be? \n";
print "Default is \"gene2go\.tab\": ";
$input_tablemaker_file = <STDIN>;
chomp ($input_tablemaker_file);

unless ($input_tablemaker_file =~ /\w+/)
{
    $input_tablemaker_file = "gene2go.tab";
}

if ($input_tablemaker_file eq "gene2go.tab") 
{
    $primary_outfile = ("wormpep2go");
    $error_list      = ("wormpep2go.errorlist");
}
else 
{
    $primary_outfile = $input_tablemaker_file . ".wormpep2go";
    $error_list      = $input_tablemaker_file . ".wormpep2go.errorlist";
}

open (INPUT_TABLE_FILE, "$input_tablemaker_file") 
    || die "Input Sanger Centre file $input_tablemaker_file not found. $!\n";
open (PRIMARY_OUTFILE, ">$primary_outfile") 
    || die "Couldn't open primary output file $primary_outfile. $!\n";
open (ERROR_LIST, ">$error_list") 
    ||  die "Couldn't open error list $error_list. $!\n";

print "\n";
print "The input file is:       $input_tablemaker_file\n";
print "    the output file is:  $primary_outfile\n";
print "    the error file is:   $error_list\n";
print "\n";

# Must process stuff like this:

# "B0261.1"	"WP:CE07704"		
# "B0261.2a"	"WP:CE07705"	"INTERPRO:IPR000403"	"GO:0004428"
# "B0261.2a"	"WP:CE07705"	"INTERPRO:IPR003151"	
# "B0261.2a"	"WP:CE07705"	"INTERPRO:IPR003152"	
# "B0261.2b"	"WP:CE28813"		

# into stuff like this:

# B0261.2a	WP:CE07705	INTERPRO:IPR000403	GO:0004428

while (<INPUT_TABLE_FILE>) 
{
    # The following conditional excludes lines without GO terms.
    if ($_ =~ /^\"(\S+)\"\t\"(\S+)\"\t\"(\S+)\"\t\"(\S+)\"/) 
    {
        $cds_literal = $1;
        $wp_id = $2;
        $interpro_code = $3;
        $go_term_literal = $4;

        # Clean up Sanger's lower-case CDS names.
        $cds_literal =~ tr/a-z/A-Z/;

        # Print all-caps CDS column. 
        print PRIMARY_OUTFILE "$cds_literal";
        print PRIMARY_OUTFILE "\t";

        # Print "WP:..." column.
        print PRIMARY_OUTFILE "$wp_id";
        print PRIMARY_OUTFILE "\t";

        # Print "INTERPRO:..." column.
        print PRIMARY_OUTFILE "$interpro_code";
        print PRIMARY_OUTFILE "\t";

        # Print "GO:..." column.
        print PRIMARY_OUTFILE "$go_term_literal";
        print PRIMARY_OUTFILE "\t";

        # Don't forget to newline separate lines!
        print PRIMARY_OUTFILE "\n";
    }
    elsif ($_ =~ /^\"(\S+)\"\t\"\S+\"\t\"\S+\"\t\"INTERPRO:IPR(\d{6})\"/) 
    {
        $cds_literal   = $1;
        $interpro_code = $2;
        $cds_literal =~ tr/a-z/A-Z/;
        chomp ($cds_literal);
        chomp ($interpro_code);
        open (INPUT_INTERPRO2GO, "$input_interpro2go") 
            || die "This program requires a readable interpro2go table! $!\n";
        while (<INPUT_INTERPRO2GO>) 
        { 
            $latest_interpro2go_line = $_;
            chomp ($latest_interpro2go_line);
            foreach ($latest_interpro2go_line) 
            {
                if ($latest_interpro2go_line =~ /$interpro_code/) 
                {
                    print ERROR_LIST "$cds_literal\t$latest_interpro2go_line\n";
                }
            }
        }
        close INPUT_INTERPRO2GO;
    }
}

close INPUT_TABLE_FILE;
close PRIMARY_OUTFILE;
close ERROR_LIST;
