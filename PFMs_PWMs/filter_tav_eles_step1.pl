#!/usr/bin/env perl

# filter_tav_eles_step1.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/31/2010.
# Purpose: first step of total kludge to get a useful subset of a very long motif prediction list.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my @privileged_sets  = ($ARGV[0],$ARGV[1]);
my $raw_list_file    = $ARGV[2];

my $serial_no = q{};
my $nt_seq    = q{};
my $ori1      = q{};
my $ori2      = q{};
my $go_term   = q{};
my $go_pval   = q{};
my $to_print  = q{};

my %priv_serial_nos = ();

foreach my $priv_file (@privileged_sets) { 
open my $PRIV, '<', $priv_file or die "Can't open $priv_file: $!";
    while (my $input = <$PRIV>) { 
        chomp $input;
        if ( $input =~ / \A \" [ACGT]+ \" \s+ (\d+) \s /xms ) { 
            $serial_no = $1;
            $priv_serial_nos{$serial_no} = 1;
        }
    }
    close $PRIV or die "Can't close filehandle to $priv_file: $!";
}

open my $RAW_LIST, '<', $raw_list_file or die "Can't open $raw_list_file: $!";
while (my $input = <$RAW_LIST>) { 
    chomp $input;

# Sample input lines:
# 
#     *... first interval is *not* simple \t!
# 2        CTGCGTCTC      381.9   557.5   0       2000    3.8e-09**       1.0e+00 3.81                                                    
# 3        AGACGCAGA      320.0   630     0       2000    1.0e+00 2.2e-08**       3.41                                                    
# 4        CGACACTCC      241.5   234     0       1500    4.2e-01 5.8e-01 2.49            positive regulation of growth   8.95e-08        

    if ( $input =~ / \A (\d+) 
                        \s+ 
                        ([ACGT]+) 
                        \t (?:[^\t]*\t){4}  
                        (\S+) 
                        \t 
                        (\S+) 
                        [^\t]* \t [^\t]* \t [^\t]* \t
                        ([^\t]*)
                        \t
                        ([^\t]*)
                        /xms ) { 
        $serial_no = $1;
        $nt_seq    = $2; 
        $ori1      = $3;
        $ori2      = $4;
        $go_term   = $5;
        $go_pval   = $6;
        $to_print  = 0;
        if ( $serial_no <= 30 ) { 
            $to_print = 1;
        }
        if (exists $priv_serial_nos{$serial_no}) { 
            $to_print = 1;
        }
        if ( $ori1 =~ /\*\*\z/xms ) { 
            if ( $ori2 =~ /\*\*\z/xms ) { 
                die "Motif no. $serial_no ($nt_seq) cannot be biased for both sense and antisense orientations simultaneously!\n";
            }
            $to_print = 1;
        }
        if ( $ori2 =~ /\*\*\z/xms ) {
            $to_print = 1;
        }
        if ( ( $go_term =~ /\S/xms ) and ( $go_pval =~ /\S/xms ) ) { 
            if (! ( looks_like_number $go_pval ) ) {
                die "Putative GO p-value $go_pval does not look numerical.\n";
            }
            $to_print = 1;
        }
        if ($to_print) { 
            print "$input\n";
        }
    }
}
close $RAW_LIST or die "Can't close filehandle to $raw_list_file: $!";

