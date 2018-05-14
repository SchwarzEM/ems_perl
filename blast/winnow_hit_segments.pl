#!/usr/bin/perl -w

# winnow_hit_segments.pl
# Erich Schwarz, 11/14/01.
# Designed to sluice out vast numbers of redundant psi-blast hits
#   and save only the roughly best ones for use in another round.
# See at the end for detailed comments.

# Datestamping is used to ensure that a given analyses' files
#   are kept sorted out; two different runs of the script will
#   automatically generate distinctly stamped files and directories.

$date = `date +"%s"`;
chomp($date);

# The list starts with biggest files; once we've seen one 
#   gi number, we don't want any more versions of it.  So
#   we set up a list of has-been reject-now gi numbers; it 
#   begins empty.

@previously_seen_sequence="";

# Set up lists, directories.

$filename="";

$discard_directory   =$date.".discard_seq_dir";
$use_directory       =$date.".seqs_to_use_dir";

$discard_list        =$date.".discard_seq_list";
$use_list            =$date.".seq_to_use_list";
$select_hit_number   =$date.".select_hit_number";

$discard_list_1      =$discard_list.".1";
$use_list_1          =$use_list.".1";
$select_hit_number_1 =$select_hit_number.".1";

system "mkdir $discard_directory";
system "mkdir $use_directory";

open DISCARD_LIST_1, ">$discard_list_1";
open USE_LIST_1, ">$use_list_1";
open HIT_LIST_1, ">$select_hit_number_1";

# Scan the working directory for files in NUMBER_aa_Res1-Res2 format.
# Mark biggest one by "ls -S" for keeping.

foreach (`ls -S | grep "_aa_"`)
{
    $filename = $_;
    if ($filename =~ /^(\d+)_aa_\d+-\d+$/) 
    {
        $discarded="no";
        for ($i = 0; $i <= $#previously_seen_sequence; $i += 1) 
        { 
            if ($1 eq $previously_seen_sequence[$i]) 
            {
            print DISCARD_LIST_1 $filename;
            $discarded="yes";
            }
        }
        if ($discarded eq "no") 
        {
            print HIT_LIST_1 "$1\n";
            push (@previously_seen_sequence, $1);
            print USE_LIST_1 $filename;
        }
    }
}

close DISCARD_LIST_1;
close USE_LIST_1;
close HIT_LIST_1;

# Tidy up the first output files.

system "sort $discard_list_1 > $discard_list";
system "sort $use_list_1 > $use_list";
system "sort $select_hit_number_1 > $select_hit_number";

system "rm $discard_list_1";
system "rm $use_list_1";
system "rm $select_hit_number_1";

# Next, move big herds of files into the rubbish, and
#   COPY (not just move) good files into a good directory.

open DISCARD_LIST, "$discard_list";
while (<DISCARD_LIST>) 
{
    $to_toss=$_;
    chomp($to_toss);
    system "mv", $to_toss, $discard_directory;
}

close DISCARD_LIST;

open USE_LIST, "$use_list";
while (<USE_LIST>) 
{ 
    $to_keep=$_;
    chomp($to_keep);
    system "cp", $to_keep, $use_directory;
}
close USE_LIST;

# for each discard list filename
# mv discard list filename into discard directory

# for each use list filename
# cp use list filename into use directory

# Given lots of files like this:

# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.extract
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.ensembl_blast.sh
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.ens_mouse_blast.sh
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.2d_30xpsi-blast.sh
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.filter.sh
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.2d_blast.sh
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.rough_filterscript
# 9961345_aa_1-385
# 9961345_aa_1-370
# 9961345_aa_1-369
# 9961345_aa_2-370
# 9961345_aa_5-370
# 9961345_aa_10-372
# 9961345_aa_9-371
# 9961345_aa_15-371
# 9961345_aa_18-365
# 9961345_aa_15-354
# 9961345_aa_6-336
# 9961345_aa_2-224
# 9961345_aa_11-231
# 9961345_aa_5-370.psi-blast-x30.nr_3.163e-04.align-6.hit_gi_nos

# pick the topmost file in the ls for each gi number:

# 9961345_aa_1-385

# and move the rest:

# 9961345_aa_1-370
# 9961345_aa_1-369
# 9961345_aa_2-370
# 9961345_aa_5-370
# 9961345_aa_10-372
# 9961345_aa_9-371
# 9961345_aa_15-371
# 9961345_aa_18-365
# 9961345_aa_15-354
# 9961345_aa_6-336
# 9961345_aa_2-224
# 9961345_aa_11-231

# into a directory for safe removal later.
