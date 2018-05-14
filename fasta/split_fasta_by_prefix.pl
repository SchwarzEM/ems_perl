#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile   = $ARGV[0];
my $nametext = $ARGV[1];

if ( (! $infile) or (! $nametext ) ) { 
    die "Format: split_fasta_by_prefix.pl [input FASTA file] [main text of output names]\n";
}

my %prefix2seen = ();

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) { 
    chomp $input;
    if ( $input =~ /\A > ([^_\.\s]+) (?:_|\.) /xms ) { 
        my $prefix = $1;
        $prefix2seen{$prefix} = 1;
    }
}
close $INFILE;

my @prefixes = sort keys %prefix2seen;

foreach my $prefix (@prefixes) { 
    my $do_print = 0;

    my $outfile = $prefix . $nametext . '.fa'; 
    $outfile    = safename($outfile);
    open my $OUTFILE, '>', $outfile;

    open  $INFILE, '<', $infile;

    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A > /xms ) { 
            if ( $input =~ /\A > $prefix _ /xms ) { 
                print $OUTFILE "$input\n";
                $do_print = 1;
            }
            else { 
                $do_print = 0;
            }
        }
        elsif ($do_print) { 
            print $OUTFILE "$input\n";
        }
    }
    close $INFILE;
    close $OUTFILE;
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


