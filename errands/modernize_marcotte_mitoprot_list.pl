#!/usr/bin/perl -w

# modernize_marcotte_mitoprot_list.pl
# Erich Schwarz, 2/1/02

# Tiny hack to update Marcotte mitochondrial protein lists.

print "Input file?: ";
$input_file = <STDIN>;
chomp ($input_file);
$output_file = $input_file . ".updated_mito";

open (INPUT_FILE, "$input_file")    || die "Can't open input file $input_file. $!\n";
open (OUTPUT_FILE, ">$output_file") || die "Can't open output file $output_file. $!\n";
while (<INPUT_FILE>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line eq "C12D12.2") 
    {
        print OUTPUT_FILE "C12D12.2A\n";
        print OUTPUT_FILE "C12D12.2B\n";
    }
    elsif ($input_line eq "C33H5.18") 
    {
        print OUTPUT_FILE "C33H5.18A\n";
        print OUTPUT_FILE "C33H5.18B\n";
    }
    elsif ($input_line eq "C53H9.2") 
    {
        print OUTPUT_FILE "C53H9.2A\n";
        print OUTPUT_FILE "C53H9.2B\n";
    }
    elsif ($input_line eq "F14F4.3") 
    {
        print OUTPUT_FILE "F14F4.3A\n";
        print OUTPUT_FILE "F14F4.3B\n";
    }
    elsif ($input_line eq "F22B7.5") 
    {
        print OUTPUT_FILE "F22B7.5A\n";
        print OUTPUT_FILE "F22B7.5B\n";
    }
    elsif ($input_line eq "F25B5.6") 
    {
        print OUTPUT_FILE "F25B5.6A\n";
        print OUTPUT_FILE "F25B5.6B\n";
    }
    elsif ($input_line eq "F26E4.13") 
    {
        print OUTPUT_FILE "F26E4.10\n";
    }
    elsif ($input_line eq "F39C12.2") 
    {
        print OUTPUT_FILE "F39C12.2A\n";
        print OUTPUT_FILE "F39C12.2B\n";
        print OUTPUT_FILE "F39C12.2C\n";
    }
    elsif ($input_line eq "F43H9.2") 
    {
        print OUTPUT_FILE "F43H9.2A\n";
        print OUTPUT_FILE "F43H9.2B\n";
    }
    elsif ($input_line eq "F52E4.1") 
    {
        print OUTPUT_FILE "F52E4.1A\n";
        print OUTPUT_FILE "F52E4.1B\n";
    }
    elsif ($input_line eq "F54H12.1") 
    {
        print OUTPUT_FILE "F54H12.1A\n";
        print OUTPUT_FILE "F54H12.1B\n";
    }
    elsif ($input_line eq "F57B9.4") 
    {
        print OUTPUT_FILE "F57B9.4A\n";
        print OUTPUT_FILE "F57B9.4B\n";
    }
    elsif ($input_line eq "M02F4.4") 
    {
        print OUTPUT_FILE "C33D12.2\n";
    }
    elsif ($input_line eq "T25D3.1") 
    {
        print OUTPUT_FILE "T25D3.3\n";
    }
    elsif ($input_line eq "Y23H5A.7") 
    {
        print OUTPUT_FILE "Y23H5A.7A\n";
        print OUTPUT_FILE "Y23H5A.7B\n";
    }
    else 
    {
        print OUTPUT_FILE "$input_line\n";
    }
}

close INPUT_FILE;
close OUTPUT_FILE;
