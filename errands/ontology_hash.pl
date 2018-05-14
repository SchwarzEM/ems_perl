#!/usr/bin/perl -w

# Program: ontology_hash.pl
#    Erich Schwarz, 12/14/01, with minor tweaks on 1/24/02.
# 
# Purpose: Create hash of "GO:---" to "C","F", and "P" for public gene 
# assignment table with one line per GO term.

print "What will input component ontology file be?\n";
print "Default is \"component.ontology\": ";
$component_infile = <STDIN>;
chomp ($component_infile);
unless ($component_infile =~ /\w+/)
{
    $component_infile = "component.ontology";
}

$component_outfile = ($component_infile . ".hash");

# 2. Open the files or die.

open (INFILE, "$component_infile") || die "Wormpep file $component_infile not found\n";
open (OUTFILE, ">$component_outfile.1") || die "Couldn't open file $component_outfile.1\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
    $line_to_be_extracted = $_;
    {
        while ($line_to_be_extracted =~ /(GO:\d{7})/) 
        {
            print OUTFILE "$1\n";
            $line_to_be_extracted =~ s/(GO:\d{7})//;
	}
    }
}

close INFILE;
close OUTFILE;

# Output file needs to be sort-ed and uniq-ed; it'd be 
#    easiest (for me at my current, rudimentary skill level) 
#    to do this by system calls.

system "sort $component_outfile.1 | uniq - > $component_outfile.2";

open (INFILE, "$component_outfile.2") || die "Wormpep file $component_outfile.2 not found\n";
open (OUTFILE, ">$component_outfile") || die "Couldn't open file $component_outfile\n";

while (<INFILE>)  
{
    print OUTFILE $_;
    print OUTFILE "C\n";
}

# Cleanup:

system "rm $component_outfile.1; rm $component_outfile.2";


print "What will input function ontology file be?\n";
print "Default is \"function.ontology\": ";
$function_infile = <STDIN>;
chomp ($function_infile);
unless ($function_infile =~ /\w+/)
{
    $function_infile = "function.ontology";
}

$function_outfile = ($function_infile . ".hash");

# 2. Open the files or die.

open (INFILE, "$function_infile") || die "Wormpep file $function_infile not found\n";
open (OUTFILE, ">$function_outfile.1") || die "Couldn't open file $function_outfile.1\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
    $line_to_be_extracted = $_;
    {
        while ($line_to_be_extracted =~ /(GO:\d{7})/) 
        {
            print OUTFILE "$1\n";
            $line_to_be_extracted =~ s/(GO:\d{7})//;
	}
    }
}

close INFILE;
close OUTFILE;

# Output file needs to be sort-ed and uniq-ed; it'd be 
#    easiest to do this by system calls.

system "sort $function_outfile.1 | uniq - > $function_outfile.2";

open (INFILE, "$function_outfile.2") || die "Wormpep file $function_outfile.2 not found\n";
open (OUTFILE, ">$function_outfile") || die "Couldn't open file $function_outfile\n";

while (<INFILE>)  
{
    print OUTFILE $_;
    print OUTFILE "F\n";
}

# Cleanup:

system "rm $function_outfile.1; rm $function_outfile.2";


print "What will input process ontology file be?\n";
print "Default is \"process.ontology\": ";
$process_infile = <STDIN>;
chomp ($process_infile);
unless ($process_infile =~ /\w+/)
{
    $process_infile = "process.ontology";
}

$process_outfile = ($process_infile . ".hash");

# 2. Open the files or die.

open (INFILE, "$process_infile") || die "Wormpep file $process_infile not found\n";
open (OUTFILE, ">$process_outfile.1") || die "Couldn't open file $process_outfile.1\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) 
{
    $line_to_be_extracted = $_;
    {
        while ($line_to_be_extracted =~ /(GO:\d{7})/) 
        {
            print OUTFILE "$1\n";
            $line_to_be_extracted =~ s/(GO:\d{7})//;
	}
    }
}

close INFILE;
close OUTFILE;

# Output file needs to be sort-ed and uniq-ed; it'd be 
#    easiest to do this by system calls.

system "sort $process_outfile.1 | uniq - > $process_outfile.2";

open (INFILE, "$process_outfile.2") || die "Wormpep file $process_outfile.2 not found\n";
open (OUTFILE, ">$process_outfile") || die "Couldn't open file $process_outfile\n";

while (<INFILE>)  
{
    print OUTFILE $_;
    print OUTFILE "P\n";
}

# Cleanup:

system "rm $process_outfile.1; rm $process_outfile.2";


# Note that to get $component/function/process_outfile back in as a hash, 
# they have be to be treated as lists.


system "touch ontology.pre-hash";
system "cat $component_outfile >> ontology.pre-hash";
system "cat $function_outfile >> ontology.pre-hash";
system "cat $process_outfile >> ontology.pre-hash";

open (ONTOLOGY_PRE_LIST, "ontology.pre-hash") || die "can't open ontology.pre-hash\n";
open (ONTOLOGY_LIST, ">ontology.hash") || die "can't open ontology.hash\n";

$ontology_term = "";

while (<ONTOLOGY_PRE_LIST>) 
    {
        $ontology_term = $_;
        chomp($ontology_term);
        $ontology_term =~ s/\:/\_/g;
        print ONTOLOGY_LIST "$ontology_term\n";
    }

close ONTOLOGY_PRE_LIST;
close ONTOLOGY_LIST;

system "rm {$component_outfile,$function_outfile,$process_outfile,ontology.pre-hash}";

print "\n";
print "The ontology hash file is \"ontology.hash\"\n";
print "\n";
