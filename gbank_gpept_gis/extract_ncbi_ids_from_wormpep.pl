#!/usr/bin/perl -w

# extract_ncbi_ids_from_wormpep.pl
# Erich Schwarz, 3/14/02.

print "Input wormpepN file: ";
$input_wp = <STDIN>;
chomp ($input_wp);
$output_list = $input_wp . ".id_numbers";
$error_list  = $input_wp . ".no_ids_list";

open (INPUT, "$input_wp") || die "Can't open $input_wp -- $!\n";
open (OUTPUT, ">$output_list") || die "Can't open $output_list -- $!\n";
open (ERRORS, ">$error_list") || die "Can't open $error_list -- $!\n";

while (<INPUT>) 
{
    $input_line = $_;
    chomp ($input_line);
    if ($input_line =~ /^>.+protein_id:(\S+)\s*.*$/) 
    {
        print OUTPUT "$1\n";
    }
    elsif ($input_line =~ /^>(.+)/)
    {
        $error_line = $1;
        chomp ($error_line);
        print ERRORS "$error_line\n";
    }
}

close INPUT;
close OUTPUT;

$output_linecount = `wc -l $output_list`;
chomp ($output_linecount);
$errors_linecount = `wc -l $error_list`;
chomp ($errors_linecount);

print "\n";
print "$output_linecount protein id. nos. were printed to: $output_list \n";
print "$errors_linecount protein id. nos. were printed to: $error_list \n";
print "\n";

# Input format of header lines:
#
# >2L52.1 CE20433    status: TR:Q9XWB3 protein_id:CAA21776.1
# >3R5.1 CE24758    status:Predicted TR:Q9XWB2 protein_id:CAB76729.1
# >4R79.2 CE19650   Ras family status:Predicted TR:Q9XXA4 protein_id:CAA20282.1
# >AC3.2 CE05132   UDP-glucuronosyltransferase status:Predicted TR:Q17399 protein_id:CAA94866.1
# >AC3.3 CE05133  locus:pqn-1  status:Predicted TR:Q17400 protein_id:CAA94867.1
