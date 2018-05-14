#!/usr/bin/perl

# process_rnai2go_tab.pl
# Erich Schwarz, 1/27/02.
# 
# To edit rnai2go.tab files by deleting "WT" lines
#    and converting "Oth" phenotypes to RNAi result names.

print "Input RNAi-to-GO table to process [default: \"rnai2go.tab\"] ";
$input_file = <STDIN>;
chomp ($input_file);
unless ($input_file =~ /\w+/) 
{
    $input_file = "rnai2go.tab";
}

$rough_input_file = $input_file . ".pre-processed";
$bad_rnai_result_file = $input_file . ".bad_rnai_result_file";

system "mv $input_file $rough_input_file";
open (ROUGH_INPUT, "$rough_input_file") || die "$rough_input_file not opened. $!\n";
open (OUTPUT, ">$input_file");
open (BAD_RNAIS, ">$bad_rnai_result_file");

# Sample input:
# 
# "JA:B0025.1"	"WT"	"B0025.1a"	"WP:CE24759"
# "JA:B0025.1"	"WT"	"B0025.1b"	"WP:CE24760"
# 
# "JA:C47B2.4"	"Oth"	"C47B2.4"	"WP:CE17564"
# "JA:F32H2.5"	"Oth"	"F32H2.5"	"WP:CE09880"
# "JA:Y71F9A_279.b"	"Oth"	"Y71F9AM.5"	"WP:CE26780"

while (<ROUGH_INPUT>)
{
    $input_line = $_;
    chomp ($input_line);
    unless ($input_line =~ /^\"(.+)\"\t\"(\w{2,3})\"(.+)/)
    {
        print BAD_RNAIS "$input_line\n";
    }
    if ($input_line =~ /^\"(.+)\"\t\"(\w{2,3})\"(.+)/) 
    {
        $rnai_result = $1;
        $rnai_phenotype = $2;
        $trailing_result = $3;
        if ($rnai_phenotype =~ /Oth/) 
        {
            $hash_friendly_rnai_result = $rnai_result;
            $hash_friendly_rnai_result =~ s/:/_/g;
            $hash_friendly_rnai_result =~ s/\./_/g;
            print OUTPUT "\"$rnai_result\"\t\"$hash_friendly_rnai_result\"$trailing_result\n";
        }
        elsif ($rnai_phenotype ne "WT")
        {
            print OUTPUT "$input_line\n";
        }
    }
}

close ROUGH_INPUT;
close OUTPUT;
