#!/usr/bin/perl
#
# check gene_association file
# version: $Revision: 1.8 $
# date: $Date: 2004/03/14 17:09:02 $
#
# specification of the gene_association format is defined at:
#   http://www.geneontology.org/GO.annotation.html#file
#
# Usage:  check-gene-association.pl [-d]
#         the -d option causes a line by line report of errors identified
#
# Maintained by the Gene Ontology Consortium
#  contact

use strict;

# process command line arguments
our $opt_d;
use Getopt::Std;

getopts('d');

# TRUE is the user wants the details report
# otherwise just the summary is provided
my $detail = defined($opt_d);

# number of line in the input file
my $linenum = 0;

# current line of text
my $line = "";

# array of errors, column number is the index
my @errors = ();

# total errors, column specific errors and line errors
my $totalerr = 0;

# errors with the whole line
my $lineerr = 0;

# Current year for date test
my $currentyear = (localtime())[5] + 1900;

# earliest time annotations were found
# this is for the date check
use constant MINYEAR => 1985;

# defined information about each column in the gene_association files
my @column = ();

#
# Column positions 
#
use constant DB => 0;
use constant DB_OBJECT_ID => 1;
use constant DB_OBJECT_SYMBOL => 2;
use constant QUALIFIER => 3;
use constant GOID => 4;
use constant REFERENCE => 5;
use constant EVIDENCE => 6;
use constant WITH => 7;
use constant ASPECT => 8;
use constant DB_OBJECT_NAME => 9;
use constant DB_OBJECT_SYNONYM => 10;
use constant DB_OBJECT_TYPE => 11;
use constant TAXON => 12;
use constant DATE => 13;
use constant ASSIGNED_BY => 14;

# Number of TAB delimited columns in file
use constant COLNUM => 15;

#
# Definition of positions in column array
#
use constant LABEL => 0;
use constant CARDINAL => 1;
use constant SPECIAL => 2;

#
# Column information:  Name, Check Cardinality equals 1, Special Check Included
#
$column[DB] = ["DB", 1, 1];
$column[DB_OBJECT_ID] = ["DB_Object_ID", 1, 0];
$column[DB_OBJECT_SYMBOL] = ["DB_Object_Symbol", 1, 0];
$column[QUALIFIER] = ["Qualifier", 0, 1];
$column[GOID] = ["GOID", 1, 1];
$column[REFERENCE] = ["DB:Reference", 1, 1];
$column[EVIDENCE] = ["Evidence", 1, 1];
$column[WITH] = ["With", 0, 0];
$column[ASPECT] = ["Aspect", 1, 0];
$column[DB_OBJECT_NAME] = ["DB_Object_Name", 0, 0];
$column[DB_OBJECT_SYNONYM] = ["DB_Object_Synonym", 0, 0];
$column[DB_OBJECT_TYPE] = ["DB_Object_Type", 1, 1];
$column[TAXON] = ["Taxon", 1, 1];
$column[DATE] = ["Date", 1, 1];
$column[ASSIGNED_BY] = ["Assigned_by", 1, 0];

#
# Evidence Codes
#
my %evicodes = ( IC => 1,  IDA => 1, IEA => 1, IEP => 1, IGI => 1, IMP => 1, 
		 IPI => 1, ISS => 1, NAS => 1, ND => 1, TAS => 1, NR => 1 );

#
# Object Types
#
my %objtypes = ( gene => 1, transcript => 1, protein => 1,
		 protein_structure => 1, complex => 1 );

#
# Database abbreviations
#
my %dbnames = (  cgen => 1, ddb => 1, ensembl => 1, fb => 1, gr => 1, genedb_gmorsitans => 1,
		 genedb_lmajor => 1, genedb_pfalciparum => 1, genedb_spombe => 1,
		 genedb_tbrucei => 1, mgi => 1, pdb => 1, rgd => 1, sgd => 1, tair => 1,
		 tigr_ath1 => 1, tigr_cmr => 1, tigr_tgi => 1, tigr_tba1 => 1, uniprot => 1,
		 vida => 1, wb => 1, zfin => 1 ); 

# begin input loop
while ( defined($line = <>) ) {
    chomp $line;
    $linenum++;

# skip comment lines
    next if ($line =~ m/^\!/);

# blank line?
    if ( $line eq "" ) {
	&checkwarn ("$linenum: BLANK line, these should be deleted are start with an \'\!\'\n");
	$lineerr++;
	$totalerr++;
	next;
    }

# split TAB delimited columns
    my @cols = split(/\t/, $line);

    if ( @cols ne COLNUM) {
	&checkwarn ("$linenum: Too few or too many columns on this line, found " . @cols . ". There should be " . COLNUM . ". Line skipped.\n");
	# increment error counters
	$lineerr++;
	$totalerr++;
	next;
    }

# loop through all the columns on this line of input
    for (my $cnum=0; $cnum < @column; $cnum++) {

# Any leading or trailing spaces?
	my $value = $cols[$cnum];
	if ( ($value =~ m/^\s/) || ($value =~ m/\s$/) ) {
	    &checkwarn ($linenum . ": " . $column[$cnum][LABEL] . " column=" . ($cnum + 1) . " leading or trailing white space: \"" . $value . "\"\n");
	    # error to have leading or trailing spaces, remove them and continue
	    $cols[$cnum] =~ s/^\s+//;
	    $cols[$cnum] =~ s/\s+$//;
	    # increment error counters
	    $errors[$cnum]++;
	    $totalerr++;
	}

# Check Cardinality
	if ($column[$cnum][CARDINAL]) {
	    if ($cols[$cnum] eq "") {
		&checkwarn ($linenum . ": " . $column[$cnum][LABEL] . " column=" . ($cnum + 1) . " cardinality should equal 1, found 0: \"" . $cols[$cnum] . "\"\n");
		$errors[$cnum]++;
		$totalerr++;
	    }
	    my @field = split(/\t/, $cols[$cnum]);
	    if ( @field > 1 ) {
		&checkwarn ($linenum . ": " . $column[$cnum][LABEL] . " column=" . ($cnum + 1) . " cardinality should equal 1, found > 1: \"" . $cols[$cnum] . "\"\n");
		$errors[$cnum]++;
		$totalerr++;
	    }
	}

#
# Specific Checks
#

	if ($column[$cnum][SPECIAL]) {

# Was a valid DB abbreviation used
	    if ($cnum == DB) {
		unless ($dbnames{ lc($cols[DB]) }) {
		    &checkwarn ("$linenum: " . $column[DB][LABEL] . " column=" . (DB + 1) . " allowed database abbreviation not correct, found \"" . $cols[DB] . "\"\n");
		    $errors[DB]++;
		    $totalerr++;
		}
	    }

# Qualifier Column on NOT and contributes_to
	    if ($cnum == QUALIFIER) {
		my @field = split(/\|/, $cols[QUALIFIER]);
		foreach my $value (@field) {
		    if ( (! $value =~ m/NOT/i) && (! $value =~ m/contributes_to/i) ) {
			&checkwarn ($linenum . ": " . $column[QUALIFIER][LABEL] . " column=" . (QUALIFIER + 1) . " only NOT and contributes_to allowed \"" . $cols[QUALIFIER] . "\"\n");
			$errors[QUALIFIER]++;
			$totalerr++;
		    }
		}
	    }

# Aspect only one of P, F or C
	    if ($cnum == ASPECT) {
		unless ( ($cols[ASPECT] eq 'P') || ($cols[ASPECT] eq 'F') || ($cols[ASPECT] eq 'C') ) {
		    &checkwarn ("$linenum: " . $column[ASPECT][LABEL] . " column=" . (ASPECT + 1) . " only P, F, or C allowed \"" . $cols[ASPECT] . "\"\n");
		    $errors[ASPECT]++;
		    $totalerr++;
		}
	    }

# Was a valid Evidence code used
	    if ($cnum == EVIDENCE) {
		unless ($evicodes{ $cols[EVIDENCE] }) {
		    &checkwarn ("$linenum: " . $column[EVIDENCE][LABEL] . " column=" . (EVIDENCE + 1) . " allowed evidence codes not present, found \"" . $cols[EVIDENCE] . "\"\n");
		    $errors[EVIDENCE]++;
		    $totalerr++;
		}
	    }

# Was a valid Object type provided
	    if ($cnum == DB_OBJECT_TYPE) {
		unless ($objtypes{ lc($cols[DB_OBJECT_TYPE]) }) {
		    &checkwarn ("$linenum: " . $column[DB_OBJECT_TYPE][LABEL] . " column=" . (DB_OBJECT_TYPE + 1) . " allowed type not present, found \"" . $cols[DB_OBJECT_TYPE] . "\"\n");
		    $errors[DB_OBJECT_TYPE]++;
		    $totalerr++;
		}
	    }

# Basic format check for reference ids
	    if ($cnum == REFERENCE) {
		my @field = split(/\t/, $cols[REFERENCE]);
		foreach my $value (@field) {
		    # WB uses [ & ] in their reference IDs
		    unless ( $value =~ m/\w+\:[\[\]\w]+/ ) {
			&checkwarn ("$linenum: " . $column[REFERENCE][LABEL] . " column=" . (REFERENCE + 1) . " format of reference not DB:REFID \"" . $cols[REFERENCE] . "\"\n");
			$errors[REFERENCE]++;
			$totalerr++;
		    }
		}
	    }

# Taxon string must start with taxon:
	    if ($cnum == TAXON) {
		my @field = split(/\|/, $cols[TAXON]);
		foreach my $value (@field) {
		    unless ( $value =~ m/^taxon:/ ) {
			&checkwarn ("$linenum: " . $column[TAXON][LABEL] . " column=" . (TAXON + 1) . " must start with taxon: \"" . $cols[TAXON] . "\"\n");
			$errors[TAXON]++;
			$totalerr++;
		    }
		}
	    }

# GOID string must start with GO:
	    if ($cnum == GOID) {
		unless ( $cols[GOID] =~ m/^GO:/ ) {
		    &checkwarn ("$linenum: " . $column[GOID][LABEL] . " column=" . (GOID + 1) . " must start with GO: \"" . $cols[GOID] . "\"\n");
		    $errors[GOID]++;
		    $totalerr++;
		}
	    }

# Check Date in proper format, YYYYMMDD
# arbitarily define the MINYEAR that makes sense
	    if ($cnum == DATE) {
		if ($cols[DATE] =~ m/(\d\d\d\d)(\d\d)(\d\d)/) {
		    if ( ($1 > $currentyear) || ($1 < MINYEAR) || ($2 > 12) || ($3 > 31) ) {
			&checkwarn ("$linenum: " . $column[DATE][LABEL] . " column=" . (DATE + 1) . " bad date format \"" . $cols[DATE] . "\"\n");
			$errors[DATE]++;
			$totalerr++;
		    }
		} elsif ($cols[DATE] ne "") {
		    # can ignore blank columns because the cardinality check would have
		    # already reported them.
			&checkwarn ("$linenum: " . $column[DATE][LABEL] . " column=" . (DATE + 1) . " bad date format \"" . $cols[DATE] . "\"\n");
			$errors[DATE]++;
			$totalerr++;
		}
	    }
	}
    }
}

# output summary of errors

# assume TAB = 8 spaces
use constant TABWIDTH => 8;

if ($totalerr > 0) {
    print "\nNUMBER of ERRORS by COLUMN\n\n";
    print "Column Name\t\tCol#\tNumber of Errors\n";
    for (my $index=0; $index < @errors; $index++) {
	if ($errors[$index] > 0) {
	    if (length($column[$index][LABEL]) < TABWIDTH) {
		$column[$index][LABEL] .= "\t";
	    }
	    if (length($column[$index][LABEL]) < (TABWIDTH * 2)) {
		$column[$index][LABEL] .= "\t";
	    }
	    print $column[$index][LABEL] . "\t" . ($index + 1) . "\t" . $errors[$index] . "\n";
	}
    }
    print "General errors\t\t-\t" . $lineerr . "\n" if ($lineerr > 0);
    print "\nTOTAL ERRORS = " . $totalerr . "\n";
    print "TOTAL ROWS in FILE = " . $linenum . "\n\n";
} else {
    print "\nCongratulations, there are no errors.\n\n";
}


# print each error if $detail equals 1
sub checkwarn
{
    my $linenum = $_[0];

    print $linenum if ($detail);
}
