#!/usr/bin/perl -w
# 
# extract_stuff_for_Keith.pl
# 12/18/01

# Sample input:
#
# WB	WP:CE24775	nhr-23		GO:0005634	PUBMED:11159333	IEA	INTERPRO:IPR001723	C		C01H6.5	protein	taxon:6239
# WB	WP:CE24775	nhr-23		GO:0006355	PUBMED:11159333	IEA	INTERPRO:IPR001723	P		C01H6.5	protein	taxon:6239
# WB	WP:CE24775	nhr-23		GO:0007275	WB:JA:C01H6.5	IMP	WB:JA:C01H6.5	P		C01H6.5	protein	taxon:6239
# WB	WP:CE24775	nhr-23		GO:0002009	WB:JA:C01H6.5	IMP	WB:JA:C01H6.5	P		C01H6.5	protein	taxon:6239
# WB	WP:CE24775	nhr-23		GO:0009653	WB:[cgc4742]:nhr-23	IMP	WB:[cgc4742]:nhr-23	P		C01H6.5	protein	taxon:6239
# WB	WP:CE24775	nhr-23		GO:0018988	WB:[cgc4742]:nhr-23	IMP	WB:[cgc4742]:nhr-23	P		C01H6.5	protein	taxon:6239

print "Input file?: ";
$input_file=<STDIN>;
chomp($input_file);
$output_file="$input_file.Keithized";
$rough_output_file="rough.$input_file.Keithized";
open (INPUT_FILE, "$input_file");
open (ROUGH_OUTPUT_FILE, ">$rough_output_file");

while (<INPUT_FILE>) 
{
    $input_line=$_;

# problem input:
# WB	WP:CE01450	gcy-1		GO:0004672	PUBMED:11159333	IEA	INTERPRO:IPR000719	F		AH6.1	protein	taxon:6239

    if ($input_line =~ /^WB.+(WP:CE\d+)\W+\w+-\w+\W+.+(GO:\d{7})\s+.+IEA\s+(INTERPRO:IPR\d{6}).+\s+(\w+\.\w+)\s+protein/) 
    {
        print ROUGH_OUTPUT_FILE "$3\t$1\t$4\t$2\n";
    }
    elsif ($input_line =~ /^WB.+(WP:CE\d+)\W+(\w+\.\w+)\W+(GO:\d{7})\D+.+IEA\s+(INTERPRO:IPR\d{6})/)
    {
        print ROUGH_OUTPUT_FILE "$4\t$1\t$2\t$3\n";
    }
    elsif ($input_line =~ /^WB.+(WP:CE\d+)\W+\w+-\w+\W+.+(GO:\d{7})\W+.+IMP\W+WB:(\S+).+\s+(\w+\.\w+)\s+protein/)
    {
        print ROUGH_OUTPUT_FILE "$3\t$1\t$4\t$2\n";
    }
    elsif ($input_line =~ /^WB.+(WP:CE\d+)\W+(\w+\.\w+)\W+(GO:\d{7})\W+.+IMP\W+WB:(\S+)\s+/) 
    {
        print ROUGH_OUTPUT_FILE "$4\t$1\t$2\t$3\n";
    }

# problem input:
# WB				GO:0007345	WB:KK:B0250.1	IMP	WB:KK:B0250.1	P		B0250.3	protein	taxon:6239

    elsif ($input_line =~ /^WB.+(GO:\d{7})\W+.+IMP\W+WB:(\S+)\s+.+\s+(\w+\.\w+)\s+protein/)
    {
        print ROUGH_OUTPUT_FILE "$2\t[what? no WP:CE number?]\t$3\t$1\n";
    }


    elsif ($input_line =~ /^WB.+/) 
    {
        print "\n";
        print "Couldn't read: \n";
        print "$input_line\n";
        print "because of stupid pattern-scan bug.\n";
        print "\n";
        die "Aiee! $!\n";
    }
}
close INPUT_FILE;
close ROUGH_OUTPUT_FILE;
system "sort $rough_output_file | uniq - > $output_file";
system "rm $rough_output_file";
