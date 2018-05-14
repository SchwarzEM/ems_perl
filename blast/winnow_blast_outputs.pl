#!/usr/bin/perl -w

# winnow_blast_outputs.pl
# Erich Schwarz, 11/21/01.
# Designed to sluice out large numbers of uninteresting Blast searches
#   and save only the ones that gave some sort of actual hit.

# Datestamping is used to ensure that a given analyses' files
#   are kept sorted out; two different runs of the script will
#   automatically generate distinctly stamped files and directories.

print "\n";
print "    Please enter a string that reliably defines\n";
print "    the Blast outputs to be manipulated in this\n";
print "    directory.  E.g. \"ensembl.genscan.fa_5e-02.align-6\".\n";
print "\n";
print "    The string is?:   ";

$defining_string=<STDIN>;
chomp($defining_string);

$date = `date +"%s"`;
chomp($date);

# Set up lists, directories.

$filename="";

$discard_directory   =$date.".discard_blast_output_dir";
$use_directory       =$date.".blast_outputs_to_use_dir";

$discard_list        =$date.".discard_blast_output_list";
$use_list            =$date.".blast_outputs_to_use_list";

$discard_list_1      =$discard_list.".1";
$use_list_1          =$use_list.".1";

system "mkdir $discard_directory";
system "mkdir $use_directory";

open DISCARD_LIST_1, ">$discard_list_1";
open USE_LIST_1, ">$use_list_1";

# Scan the working directory for files in NUMBER_aa_Res1-Res2 format.
# Mark biggest one by "ls -S" for keeping.

$discard="maybe";

foreach (`ls | grep $defining_string`)
{
    if ($discard eq "no") 
    {
        print USE_LIST_1 "$filename\n";
        system "cp", $filename, $use_directory;
    }
    $discard="no";
    $filename = $_;
    chomp($filename);
    open FILENAME, "$filename";
    while (<FILENAME>) 
    {
        if ($_ =~ /No hits found/)
        {
            $discard="yes";
            print DISCARD_LIST_1 "$filename\n";
            system "mv", $filename, $discard_directory;
        }
    }
    close FILENAME;
}

close DISCARD_LIST_1;
close USE_LIST_1;

# Tidy up the first output files.

system "sort $discard_list_1 > $discard_list";
system "sort $use_list_1 > $use_list";

system "rm $discard_list_1";
system "rm $use_list_1";

print "\n";
print "Job done.\n";
print "\n";
