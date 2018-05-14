#!/usr/bin/perl -w

chomp(my @files = `ls`);

foreach my $file (@files) 
{
    my $newfile = $file;
    $newfile =~ s/ /_/g;
    if (-e $newfile) 
    {
        warn "Can't rename $file to $newfile, because $newfile already exists.\n";
    }
    elsif (rename $file, $newfile)
    {  
        # success means no other action needed
    }
    else
    {
        warn "Renaming $file to $newfile failed: $!\n";
    }
}
