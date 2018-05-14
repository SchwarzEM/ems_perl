#!/usr/bin/env perl

# bed_vs_fimo.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/11/2010.
# Purpose: given a *.bed file, and a fimo hitlist of motifs against that *.bed's FASTA, print either positive or negative subset BED.

use strict;
use warnings;
use Getopt::Long;

my $positive;
my $negative;
my $help;

my $bed_file      = q{};
my $fimo_file     = q{};
my $ele_id        = q{};
my %ele_exists    = ();
my %ele_has_motif = ();

GetOptions ( 'bed=s'    => \$bed_file,
             'fimo=s'   => \$fimo_file,
             'positive' => \$positive,
             'negative' => \$negative,
             'help'     => \$help,           );

if ( $help 
     or ( (! $positive ) and (! $negative )              )
     or ( $positive and $negative                        )
     or (! (-r $bed_file )                               ) 
     or (! ( (-r $fimo_file ) or ( $fimo_file eq '-' ) ) ) ) { 
    die "Format: bed_vs_fimo.pl",
        " --bed|-b [BED file]",
        " --fimo|-f [fimo.txt] file",
        " --positive|-p [select elements with fimo hits]",
        " <or>",
        " --negative|-n [select elements *without* fimo hits]",
        "\n",
        ;
}

# Read BED for the first time, just to get IDs for elements.
open my $BED1, '<', $bed_file or die "Can't open BED file $bed_file: $!";
while (my $input = <$BED1>) { 
    chomp $input;

    # Silently ignore lines starting with 'track name="', but require parsing of all others.
    if (      ( $input !~ /\A track \s+ name = \" /xms   )
          and ( $input =~ / \A (?: \S+ \t){3} (\S+) /xms ) ) { 
        $ele_id = $1;
        $ele_exists{$ele_id} = 1;
    }

    # Die loudly if parsing fails.
    elsif ( $input !~ /\A track \s+ name = \" /xms ) { 
        die "In BED file $bed_file, can't parse input line: $input\n";
    }
}
close $BED1 or die "Can't close handle to BED file $bed_file: $!";

# Read fimo.txt file to get each subset of elements:
open my $FIMO, '<', $fimo_file or die "Can't open FIMO text file $fimo_file: $!";
while (my $input = <$FIMO>) { 
    chomp $input;

    # Silently ignore lines starting with '#', but require parsing of all others.
    if (     ( $input !~ /\A \# /xms                            ) 
         and ( $input =~ /\A \S+ \t [^\s]+_([^\s_]+) \t /xms ) ) { 
        $ele_id = $1;
        if (! ( exists $ele_exists{$ele_id} ) ) { 
            die "In FIMO file $fimo_file, can't identify element \"$ele_id\" in input line: $input\n";
        }

        # Record which elements have one or more motifs, for later BED reading and subset-producing.
        $ele_has_motif{$ele_id} = 1;
    }

    # Die loudly if parsing fails.
    elsif ( $input !~ /\A \# /xms ) { 
        die "In FIMO file $fimo_file, can't parse input line: $input\n";
    }
}
close $FIMO or die "Can't close handle to FIMO text file $fimo_file: $!";

# Read BED for second time; this time, print out results.
open my $BED2, '<', $bed_file or die "Can't open BED file $bed_file: $!";
while (my $input = <$BED2>) {
    chomp $input;

    my $head_name = q{};
    my $head_desc = q{};
    my $head_text = q{};

    my $blurb = q{};
    if ($positive) {
        $blurb = 'with_motifs';
    }
    if ($negative) {
        $blurb = 'no_motifs';
    }

    # Header line:
    if ( $input =~ / \A track \s+ name = \" ([^\"]+) \" \s+ description = \" ([^\"]+) \" (.*) \z /xms ) {
        $head_name = $1;
        $head_desc = $2;
        $head_text = $3;
        print qq{track name=\"$head_name},
               q{_},
              qq{$blurb\" description=\"$head_desc},
               q{_},
              qq{$blurb\"$head_text\n},
               ;
    }

    # Otherwise, require parsing/filtering of lines or death.
    if (      ( $input !~ /\A track \s+ name = \" /xms   )
          and ( $input =~ / \A (?: \S+ \t){3} (\S+) /xms ) ) {
        $ele_id = $1;
        if (! ( exists $ele_exists{$ele_id} ) ) {
            die "In BED file $bed_file, can't identify element \"$ele_id\" in input line: $input\n";
        } 
        if ( $positive and ( exists $ele_has_motif{$ele_id} ) ) { 
            print "$input\n";
        }
        if ( $negative and (! exists $ele_has_motif{$ele_id} ) ) {
            print "$input\n";
        }
    }
    elsif ( $input !~ /\A track \s+ name = \" /xms ) {
        die "In BED file $bed_file, can't parse input line: $input\n";
    }
}
close $BED2 or die "Can't close handle to BED file $bed_file: $!";

