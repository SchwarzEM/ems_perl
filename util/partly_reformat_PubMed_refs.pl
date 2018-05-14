#!/usr/bin/env perl

# partly_reformat_PubMed_refs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/8/2010.
# Purpose: convert (e.g) "Lin JJ, Lobo NF, Lopez JR, Malek JA ... 2008" to "Lin, J.J, Lobo, N.F., Lopez, J.R., Malek, J.A. 2008.", etc; TODO -- figure out good way to expand page-numbers.

use strict;
use warnings;

my $data_ref;
my $i = 0;

my $input_text     = q{};
my $remaining_text = q{};
my $output_text    = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \s* \z/xms ) { 
        if ($input_text) { 
            $i++;
            $data_ref->{$i} = $input_text;
            $input_text = q{};
        }
    }
    if ( $input =~ / \S /xms ) {
        $input_text .=  q{ } . $input;
    }
}

if ($input_text) {
    $i++;
    $data_ref->{$i} = $input_text;
    $input_text = q{};
}

if ( $i > 0 ) { 
    print "\n";
}

foreach my $j (1..$i) { 
    $input_text = $data_ref->{$j};
    my $part1    = q{};
    my $initials = q{};
    my $part3    = q{};
    if ( $input_text =~ /\A ([^\.]+) [ ] ([A-Z]+) (\. .*) \z /xms ) { 
        $part1    = $1;
        $initials = $2;
        $part3    = $3;
        my @init_letts = split //, $initials;
        $initials = join '.', @init_letts;
        $input_text = $part1 . ', ' . $initials . $part3;
    }
    while ( $input_text =~ /\A ([^,]+) \s ([A-Z]+) \, (.+) /xms ) {
        my ($surname, $initials);
        $surname        = $1;
        $initials       = $2;
        $remaining_text = $3;
        $output_text .= "$surname, ";
        while ( $initials =~ / ([A-Z]{1}) (.*) /xms ) {
            my $front_text = $1;
            $initials = $2;
            $output_text .= "$front_text.";
            $initials = $2;
        }
        $output_text .= q{,};
        $input_text = $remaining_text;
    }
    $remaining_text =~ s/[ ]+/ /g;
    $remaining_text =~ s/ ,/,/g;

    # Spieth, J. WormBase: new content and better access. Nucleic Acids Res. 2007 Jan;35(Database issue):D506-10.
    if ( $remaining_text =~ / \A 
                              (.+? \.) \s
                              (.+ \. \s) 
                              ( (?: 19|20) \d{2} ) 
                              \s (?: (?: Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) (?: \s \d{1,2}){0,1} ; ) 
                              ( \w+ )
                              \( .* \) 
                              : 
                              ( [a-zA-Z\d\-]+ \. )
                              \s* 
                              \z /xms ) { 
        my $last_name   = $1;
        my $front_text  = $2;
        my $year        = $3;
        my $back_text1  = $4;
        my $back_text2  = $5;
        $remaining_text = "$last_name $year. " . $front_text . q{ } . $back_text1 . q{, } . $back_text2;
        my $commas_in_leading_names = ( $output_text =~ s/,/,/g );
        if ( $commas_in_leading_names >= 4 ) { 
            $remaining_text = ' and ' . $remaining_text;
        }
        $remaining_text =~ s/[ ]+/ /g;

    }

    $output_text .= $remaining_text;

    $output_text =~ s/\A\s+//;
    $output_text =~ s/\s+\z//;
    $output_text =~ s/[ ]+/ /g;
    $output_text =~ s/ ,/,/g;

    # Export as one line, but spaced to make cut-and-paste easy.
    print "$output_text\n\n";

    # Clear text from this variable, so that it doesn't get repeated:
    $output_text = q{};
}

