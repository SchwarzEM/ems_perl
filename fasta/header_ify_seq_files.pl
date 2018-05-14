#! /usr/bin/perl -w

# header_ify_seq_files.pl
# Erich Schwarz <emsch@its.caltech.edu>, 10/07/03

# Purpose: take a bunch of non-headered plaintext sequence files and give them FASTA headers.

print "Directory of non-headered plaintext sequence files?: ";
chomp ($input_directory = <STDIN>);

chomp($date1 = `date -I`);
chomp($date2 = `date +%s`);
$date = "$date1" . "_" . "$date2";

$output_file = $input_directory . $date . ".output";

chomp(@input_file_list = `ls $input_directory`);

chdir ($input_directory) || die "Can't change to directory $input_directory\n";
foreach $input_file (@input_file_list) 
{
                                         #  unfortunately, confusing visual syntax:
    open (OUTPUT, ">>$output_file");     #  '>>' append to a file
    print OUTPUT ">$input_file\n";       #  '>'  type '>' at start of a header line
    close OUTPUT;
    system "cat $input_file >> $output_file";
}

system "mv $output_file ..";
chdir ("..") || die "Can't return to original directory\n";

