#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @input_files = @ARGV;
my $header      = "File\tChastity_status\tSample_header_line";
my $INFILE;

print "$header\n" if @input_files;

LOOP: foreach my $input_file (@input_files) {
    open $INFILE, '<', $input_file;
    if (my $input = <$INFILE>) {
        if ( $input =~ /\A [@] \S+ \s+ 1:([NY]):0:[ACGT]+ \s* \z/xms ) { 
            close $INFILE;
            open my $INFILE2, '<', $input_file;
            my $i = 0;
            my $j = 0;
            EXTENDED_TEST: while (my $input = <$INFILE2>) {
                $i++;
                $j = ($i % 4);
                if ( ( $j == 1 ) and ( $input =~ /\A [@] \S+ \s+ 1:([NY]):0:[ACGT]+ \s* \z/xms ) ) {
                     my $chastity_flag = $1;
                     if ( $chastity_flag eq 'Y' ) {
                         my $output = "$input_file\tNEEDS_CHASTITY\t$input";
                         chomp $output;
                         print "$output\n";
                         close $INFILE2;
                         last EXTENDED_TEST;
                     }
                }
            }
            next LOOP;
        }
        elsif ( $input =~ /\A [@] \S+ .* \z/xms ) { 
            my $output = "$input_file\tno_chastity_needed\t$input";
            chomp $output;
            print "$output\n";
            close $INFILE;
            next LOOP;
        }
        else {
            my $output = "$input_file\tMALFORMATTED\t$input";  
            chomp $output;
            print "$output\n";
            close $INFILE;
            next LOOP;
        }
    }
    else { 
        die "Can't read first line of $input_file\n";
    }
}

