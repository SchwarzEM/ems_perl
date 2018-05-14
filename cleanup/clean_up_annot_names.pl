#!/usr/bin/perl -w 

# clean_up_annot_names.pl
# Erich Schwarz <emsch@its.caltech.edu>, 4/2/03

# Purpose: fix or warn about names in funct. annot.s not in 'WBPerson\d+' format.

print "Annotation file whose names are to be cleaned or warned about? ";

chomp ($input_file         = <STDIN>);
$renamed_annot_file        = $input_file . ".renamed_annot_file";
$rough_output_file         = $input_file . ".rough_output_file";
$list_of_nonstandard_names = $input_file . ".list_of_nonstandard_names";

%name_table = ("Bastiani CA "         => "WBPerson48",  
               "Bastiani CA"          => "WBPerson48",  
               "Carol Bastiani"       => "WBPerson48",  
               "Chan J"               => "WBPerson1823",
               "Gobel V"              => "WBPerson204",   
               "Hodgkin JA"           => "WBPerson261",
               "Kishore R"            => "WBPerson324",
               "Kostrouchova M"       => "WBPerson344", 
               "Kramer JM"            => "WBPerson345", 
               "Lee RYN"              => "WBPerson363",  
               "Mole SE"              => "WBPerson1832",
               "Muller B"             => "WBPerson1874",  # Berndt Mueller <b.mueller@abdn.ac.uk> -- helped annotate cdl-1
               "Petcherski AG"        => "WBPerson480", 
               "Schwarz EM"           => "WBPerson567", 
               "Sternberg PW"         => "WBPerson625",
               "Kimberly Van Auken"   => "WBPerson1843");

open (INPUT, "$input_file") || die;
open (RENAMED_ANNOTS, ">$renamed_annot_file") || die;

while (<INPUT>)
{
    chomp ($input_line = $_);
    if ($input_line =~ /(.+Person_evidence \")([^"]+)(\".*)$/) 
    {
        $input_line_part_1 = $1;
        $input_line_part_2 = $2;
        $input_line_part_3 = $3;
        if (exists $name_table{"$input_line_part_2"})
        {
            $input_line_part_2 = $name_table{"$input_line_part_2"};
        }
        $input_line = $input_line_part_1 . $input_line_part_2 . $input_line_part_3;
    }
    print RENAMED_ANNOTS "$input_line\n";
}

close INPUT; close RENAMED_ANNOTS;

open (INPUT, "$renamed_annot_file") || die;
open (ROUGH_OUTPUT, ">$rough_output_file") || die;

while (<INPUT>) 
{
    chomp ($input_line = $_);
    if ($input_line =~ /.+Person_evidence \"([^"]+)\"/) 
    {
        $input_name = $1;
        unless ($input_name =~ /WBPerson\d+/) 
        {
            print ROUGH_OUTPUT "\"$input_name\"\n";
        }
    }
}

system "sort $rough_output_file | uniq - > $list_of_nonstandard_names";
system "rm $rough_output_file";
