#!/usr/bin/perl -w

# clean_up_.ace-readable.pl
# Erich Schwarz, 4/3/02

print "Input .ace-readable file: ";
$proto_ace_input = <STDIN>;
chomp ($proto_ace_input);

$ace_output = $proto_ace_input . ".ace";

open (PROTO_ACE_INPUT, "$proto_ace_input") || die "Can't open $proto_ace_input. $!\n";
open (ACE_OUTPUT, ">$ace_output") || die "Can't open $ace_output. $!\n";

$count_of_quotes = 0;
$new_count_of_quotes = 0;
while (<PROTO_ACE_INPUT>) 
{
    $input_line = $_;
    unless ($input_line =~ /^\/\//) 
    {
        $new_count_of_quotes = tr/\"//;
        $count_of_quotes = ($count_of_quotes + $new_count_of_quotes);
        chomp ($input_line);
        if (($count_of_quotes == 0) || ($count_of_quotes == 2)) 
        {
            print ACE_OUTPUT "$input_line\n";
            $count_of_quotes = 0;
            $new_count_of_quotes = 0;
        }
        elsif ($count_of_quotes == 1) 
        {
        if ($input_line =~ /\s$/) 
            {
             print ACE_OUTPUT "$input_line";
            }
        else 
            {
            print ACE_OUTPUT "$input_line";
            print ACE_OUTPUT " ";
            }
        }
        else 
        {
        print "Miscount of quotes at this line:\n";
        print "$input_line\n";
        die;
        }
    }
}

close PROTO_ACE_INPUT;
close ACE_OUTPUT;
