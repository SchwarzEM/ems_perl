#!/usr/bin/perl -w

# dot-ace-ify_batch_set.pl
# Erich Schwarz, 5/28/02

# Purpose: automatically generate a old-style .ace file
# that can be used for functional annotation of 
# a batch of genes sharing a common trait described 
# in a common text.

print "Input list of common sequences? ";
$input_sequences = <STDIN>;
chomp ($input_sequences);

# Input sequences should be sort-ed and uniq-ed, in this format:
# ...
# H15N14.1A
# H15N14.1B
# K01B6.1
# K08E3.8
# locus:abf-6
# locus:apr-1
# ...

print "Input file of evidence hashes? ";
$input_evidence = <STDIN>;
chomp ($input_evidence);

print "Input file of common text description? ";
$input_description = <STDIN>;
chomp ($input_description);

$output_file = $input_sequences . ".old-style.ace";

open (INPUT_SEQUENCES, "$input_sequences")  || die "Bogosity! $input_sequences unopenable. $!\n";
open (OUTPUT, ">$output_file")              || die "Bogosity! $output_file unopenable. $!\n";

while (<INPUT_SEQUENCES>) 
{
    $input_line = $_;

    if ($input_line =~ /locus:(.+)\s*\n/) 
    {
        $target = $1;
        $gene_name = $target;
        print OUTPUT "\n";
        print OUTPUT "Locus : \"$target\"\n";
    }

    elsif ($input_line =~ /^(\w+.\w+)([A-Z]{1})\s*\n/)
    {
        $header = $1;
        $suffix = $2;
        $gene_name = $header;
        $suffix =~ tr/A-Z/a-z/;
        $target = $header . $suffix;
        print OUTPUT "\n";
        print OUTPUT "Sequence : \"$target\"\n";
    }

    elsif ($input_line =~ /^(\w+.\w+)\s*\n/)
    {
        $target = $1;
        $gene_name = $target;
        print OUTPUT "\n";
        print OUTPUT "Sequence : \"$target\"\n";
    }

    open (EVIDENCE, "$input_evidence") || die "Bogosity! $input_evidence unopenable. $!\n";

    while (<EVIDENCE>) 
    {
        $evidence_line = $_;
        print OUTPUT "Provisional_description \"$target\" $evidence_line";
    }
    print OUTPUT "\n";
    print OUTPUT "LongText : \"$target\"\n";
    print OUTPUT "\n";

    open (DESCRIPTION, "$input_description");
    while (<DESCRIPTION>) 
    {
        $description_line = $_;
        print OUTPUT "$description_line";
    }
    close DESCRIPTION;
    print OUTPUT "\n";
    print OUTPUT "\*\*\*LongTextEnd\*\*\*\n";
    print OUTPUT "\n";
    print "The $gene_name gene was automatically annotated.\n";
}

close INPUT_SEQUENCES;
close OUTPUT;
