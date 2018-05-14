#!/usr/bin/env perl

# Program: tfa_split.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/1/2004 + 9/26/2008.
# Purpose: split .tfa files with multiple entries into many files with single entries.

use strict;
use warnings;
use Getopt::Long;

my $date    = join('.', &get_local_date());
my $outfile = q{};
my $suffix  = q{};
my $OUTFILE;

GetOptions ("suffix=s" => \$suffix); 

while (my $input = <>) {
    if ($input =~ /\A > (\S+) /xms ) {
        my $outfile_new = $1;
        if ($suffix =~ /\A \S+ \z/xms) { 
            $outfile_new = $outfile_new . '.' . $suffix;
        }
        if ($OUTFILE) { 
            close $OUTFILE or die "Can't close filehandle to $outfile: $!";
        }
        $outfile = failsafe_name($outfile_new);
        open $OUTFILE, '>', $outfile or die "Can't open outfile $outfile\n";
        print {$OUTFILE} $input;
    } 
    else {
        print {$OUTFILE} $input; 
    }
}
close $OUTFILE or die "Can't close filehandle to $outfile: $!";

# Revised, 9/18/2010.
# The earlier version of this would truncate pre-existing "\.\d+" in name stems.

sub failsafe_name {
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

sub get_local_date {
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;
}

