#!/usr/bin/env perl

# fastq2fa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/6/2010.
# Purpose: convert FASTQ files with no information in the quality-data header line ("\A\+\s*\z") to FASTA format.

use strict;
use warnings;
use Getopt::Long;

my $should_print        = 0;
my $num_seq_lines_read  = 0;
my $expect_qual_lines   = 0;
my $num_qual_lines_read = 0;

my $current_seqname   = q{};
my $new_seqname       = q{};
my $new_header        = q{};
my %observed_seqnames = ();

my $allow_repeatnames;
my $help;

GetOptions ( "allow_repeatnames" => \$allow_repeatnames,
             "help"              => \$help, );

if ($help) { 
    die "fastq2fa.pl -h|--help -a|--allow_repeatnames <input FASTQ stream>\n";
}

while (my $input = <>) { 
    chomp $input;

    # Long if-elsif[xN] series.

    if ( $should_print and $expect_qual_lines ) { 
        die "Script is badly confused: expects both printable sequence *and* quality lines!\n";
    }

    # Accept exactly as many quality lines as there have been sequence lines.
    # (But make no attempt to check that each line set had equal lengths.)
    elsif ( $expect_qual_lines ) { 
        if ( $input =~ / \A \S+ /xms ) { 
            $num_qual_lines_read++;
            if ( $num_qual_lines_read == $num_seq_lines_read ) { 
                $expect_qual_lines = 0;
            }
            if ( $num_qual_lines_read > $num_seq_lines_read ) { 
                die "Somehow overcounted qual. lines at: $input\n";
            }
        }
    }

    # If expecting printable sequence, a line starting with '+' should be a seq.-qual. header.
    # Note that this assumes that there are *no* characters after the first '+'.
    elsif ( ( $input =~ / \A \+ \s* \z /xms ) and $should_print ) { 
        $should_print      = 0;
        $expect_qual_lines = 1;
    }

    # This should be plain DNA sequence.
    elsif ( $input =~ / \A [acgtnACGTN]+ \z /xms ) { 
        if (! $should_print ) { 
            die "Can't print sequence-like line: $input\n";
        }
        if ($should_print) { 
            print "$input\n";
            # Count no. of seq. lines to later keep quality lines equal.
            $num_seq_lines_read++;
        }
    }

    # This should be a fastq header of DNA sequence.
    # Originally this selected only '([\w\.]+)' as $2 -> $new_seqname,
    #     but it had to be broadened to '([\w\-\:\.#\/]+)' to cope with wacky CIT headers.
    elsif ( ( $input =~ / \A \@ (([\w\-\:\.#\/]+).*) \z /xms ) 
            and (! $should_print ) 
            and (! $expect_qual_lines ) ) { 
        $new_header  = $1;
        $new_seqname = $2;
        if ( ( $observed_seqnames{$new_seqname} ) and (! $allow_repeatnames) ) {
            die "Two different sequences named",
                " \"$new_seqname\", with second one",
                " in input: $input\n",
                ;
        }
        $current_seqname                     = $new_seqname;
        $observed_seqnames{$current_seqname} = 1;
        $should_print                        = 1;
        $num_seq_lines_read                  = 0; 
        $num_qual_lines_read                 = 0;
        print ">$new_header\n";
    }

    # If at least one the above conditions hasn't been met, try to loudly die.
    elsif ( ( $input =~ /\S+/xms ) and $should_print ) { 
        die "Shouldn't print input line of unclear status: $input\n";
    }
}

