#!/usr/bin/perl -w 

# copy_file_list.pl
# Erich Schwarz, 11/15/01.
# Wee scriptoid to copy files to a new directory.

print "\n";
print "List of files to copy? ";
$filelist=<STDIN>;
chomp($filelist);

print "\n";
print "Directory to copy them to? ";
$destination=<STDIN>;
chomp($destination);

unless (-e $filelist && -r $filelist) 
{ 
    print "\n";
    die "The list $filelist isn't available. $!\n";
    print "\n";
}

unless (-e $destination && -w $destination)
{
    print "\n";
    die "The directory $destination isn't available. $!\n";
    print "\n";
}

foreach $filename (`cat $filelist`) 
{
    chomp($filename);
    system "cp", $filename, $destination;
}

print "\n";
print "Files listed in $filelist have been copied to $destination.\n";
print "\n";
