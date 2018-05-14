#!/usr/bin/env perl

# cds2full_name.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/23/2010.
# Purpose: given table of WB \t cds \t CGC names, and file/stream of lines starting with sequence names, print the same but starting with full WB|cds|CGC names. 

use strict;
use warnings;
use Getopt::Long;

my $infile     = q{};
my $name_table = q{};
my $help;

my $wbg     = q{};
my $cds     = q{};
my $cgc     = q{};
my %cds2wb  = ();
my %cds2cgc = ();

my $name    = q{};
my $text    = q{};

GetOptions ( 'infile=s' => \$infile,
             'names=s'  => \$name_table,
             'help'     => \$help, );

if ( $help or (! $infile ) or (! $name_table ) ) { 
    die "Format: cds2full_name.pl --infile|-i <input stream/files>  --names|-n [tab-delimited WBGene / CDS / CGG name file] --help|-h\n";
}

open my $NAMES, '<', $name_table or die "Can't open name table file $name_table: $!";
while (my $input = <$NAMES>) { 
    chomp $input;
    if ( $input =~ / \A (WBGene\d+) \t ([^\t]+) \t ([^\t]+) \s* \z /xms ) { 
        $wbg = $1;
        $cds = $2;
        $cgc = $3;
        if ( ( $cds =~ /\-/xms ) or ( $cgc !~ /\-/xms ) ) { 
            die "Apparently misformatted input line: $input\n";
        }
        $cds2wb{$cds}  = $wbg;
        $cds2cgc{$cds} = $cgc;
    }
    elsif ( $input =~ / \A (WBGene\d+) \t ([^\t]+) \s* \z /xms ) {
        $wbg = $1;   
        $cds = $2;
        if ( $cds =~ /\-/xms ) { 
            die "Apparently misformatted input line: $input\n";
        }
        $cds2wb{$cds}  = $wbg;
    }
    else { 
        die "Apparently misformatted input line: $input\n";
    }
}
close $NAMES or die "Can't close filehandle to table file $name_table: $!";

# Accept either a stream from '-' or a standard file.
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
    if ( $input =~ /\A (\S+) (\s+ .+) \z /xms ) { 
        $name = $1;
        $text = $2;
        my $outname = q{};
        if ( exists $cds2wb{$name} ) {
            $outname .= "$cds2wb{$name}|";
        }
        $outname .= "$name";
        if ( exists $cds2cgc{$name} ) { 
            $outname .= "|$cds2cgc{$name}";
        }
        print $outname, $text, "\n";
    }
    else { 
        warn "Could not process: $input\n";
        print "$input\n";
    }
}

