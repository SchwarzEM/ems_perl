#!/usr/bin/env perl

# renumber_wb_motifs.pl -- Erich Schwarz <emsch@caltech.edu>, 10/9/2012
# Purpose: given a file or stream of WBMotif text in .ace format, and a designated starting number N, renumber the motifs starting from N (default is N=1).

use strict;
use warnings;
use Getopt::Long;

my @infiles = ();
my $number  = 1;
my $comment;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'number:i'     => \$number,
             'comment'      => \$comment,
             'help'         => \$help, );

if ( (! @infiles) or (! $number) or $help ) { 
    die "Format: renumber_wb_motifs.pl\n",
        "    --infile|-i   <input stream/files>\n",
        "    --number|-n   [starting integer for renumbered output WBPmat series]\n",
        "    --comment|-c  [append old name as a comment after new numbered name]\n",
        "    --help|-h     [print this message]\n",
        ;
}


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
        if ( $input =~ / \A (Position_Matrix \s+ : \s+ \") ([^\"\s]+) (\" .*) \z /xms ) { 
            my $pre_old_name  = $1;
            my $old_name      = $2;
            my $post_old_name = $3;
            my $new_id        = sprintf "%08i", $number;
            $new_id           = 'WBPmat' . $new_id;
            my $output        = $pre_old_name . $new_id . $post_old_name;
            if ($comment) { 
                $output = $output . q{ // formerly called: } . $old_name;
            }
            print "$output\n";
            $number++;
        }
        elsif ( $input =~ / \A Position_Matrix /xms ) { 
            die "Can't parse input: $input\n";
        }
        else {
            print "$input\n";
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

