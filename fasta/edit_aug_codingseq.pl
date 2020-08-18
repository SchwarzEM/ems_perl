#!/usr/bin/env perl

# edit_aug_codingseq.pl -- Erich Schwarz <ems394@cornell.edu>, first coded on 11/16/2012; revised on 8/20/2016 to enforce exact duplication in each name.
# Purpose: correct a buggy echo.echo output of sequence names for *.aug.codingseq products of getAnnoFasta.pl (first observed in Nov. 2012; still there in Aug. 2016).

use strict;
use warnings;
use Getopt::Long;

my @infiles      = ();
my $name_pattern = q{};
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'name=s'       => \$name_pattern,
             'help'         => \$help,   );

$name_pattern ||= q{\S+?};

if ( $help or (! @infiles) or (! $name_pattern) ) { 
    die "Format: edit_aug_codingseq.pl\n",
        "    --infile|-i   <input stream/files>\n",
        "    --name|-n     [name pattern recurring in header, e.g., ", q{"Acey_2012\.08\.05_\d+"}, "; default is ", q{"\S+?"}, "]\n",
        "    --help|-h     [print this message]\n",
        ;
}

# Do special Perl quoting to make the pattern work better in a regex:
$name_pattern = qr/$name_pattern/;

foreach my $infile (@infiles) { 
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( $input =~ /\A > /xms ) { 
            if ( $input =~ /\A > ($name_pattern) \. (($name_pattern) \. g \d+ \. t \d+ \b .*) \z/xms ) { 
                my $first_tandem_repeat  = $1;
                my $revised_header       = $2;
                my $second_tandem_repeat = $3;

                if ( $first_tandem_repeat ne $second_tandem_repeat ) {
                    die "Can't parse header line, which fails to show exact identity of \"$first_tandem_repeat\" ",
                        "versus \"$second_tandem_repeat\": $input\n",
                        ;
                }

                $revised_header = '>' . $revised_header;
                print "$revised_header\n";
            }
            else { 
                die "Can't parse header line, which fails to show duplication of pattern \"$name_pattern\": $input\n";
            }
        }
        else { 
            if ( $input =~ /\S/xms ) { 
                print "$input\n";
            }
        }
    }
}

