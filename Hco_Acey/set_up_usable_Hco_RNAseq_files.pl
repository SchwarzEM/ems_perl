#!/usr/bin/env perl

# set_up_usable_Hco_RNAseq_files.pl -- Erich Schwarz <emsch@caltech.edu>, 8/20/2011.
# Purpose: given an input text from March 2011, make a reliable script and text for rebuilding RNA-seq data.

use strict;
use warnings;

my $library_name  = q{};
my $flowcell_name = q{};
my $lane_end_nums = q{};
my $lane_num      = q{};
my $end_num       = q{};

my %labels = ( lib_11344 => 'Egg',
               lib_11345 => 'L1',
               lib_11419 => 'L2-1',
               lib_11439 => 'L2-2',
               lib_11346 => 'L3_untrusted',
               lib_11441 => 'L3-1',
               lib_11347 => 'L4_female',
               lib_11348 => 'L4_male',
               lib_11445 => 'Adult_female',
               lib_11443 => 'Adult_male',
               lib_10868 => 'Prob_mixed_adults', );

my $prefix1A = '/woldlab/loxcyc/data00/solexa-sequence/flowcells/';
my $prefix1B = '/C1-152/';

# Aimed at something like:
#    /woldlab/loxcyc/data00/solexa-sequence/flowcells/42FC3AAXX/C1-152/s_1_1_sequence.txt.bz2

my $directory1 = '/home/schwarz/public_html/Haecon_2011/raw/RNAseq';

# Aimed at something like:
#    /home/schwarz/public_html/Haecon_2011/raw/RNAseq/lib_10868_42FC3AAXX.1_1_sequence.txt.bz2

my $prefix3 = 'http://woldlab.caltech.edu/~schwarz/Haecon_2011/raw/RNAseq/';

# Aimed at something like:
#    http://woldlab.caltech.edu/~schwarz/Haecon_2011/raw/RNAseq/lib_10868_42FC3AAXX.1_1_sequence.txt.bz2 

while (my $input = <>) { 

# Typical input lines:
#     http://woldlab.caltech.edu/~schwarz/Haecon_2011/raw/RNAseq/lib_10868/42FC3AAXX.1.orig/s_1_1_sequence.txt.bz2
#     http://woldlab.caltech.edu/~schwarz/Haecon_2011/raw/RNAseq/lib_10868/42FC3AAXX.1.orig/s_1_2_sequence.txt.bz2
# or
#     http://woldlab.caltech.edu/~schwarz/Haecon_2011/raw/RNAseq/lib_10868/311B7AAXX.8.orig/s_8_sequence.txt.bz2

    chomp $input;
    if ( ( $input =~ /\A \s* [#] /xms ) or ( $input =~ /\A \s* \z/xms ) ) { 
        print "$input\n";
    }

    elsif ( $input =~ / \A 
                        \s*

                        http:\/\/
                        woldlab.caltech.edu\/
                        [~]schwarz\/Haecon_2011\/raw\/RNAseq\/
                        (lib_\d+)                                       # library name
                        \/
                        ([A-Z0-9]+)                                     # flowcell name
                        \. \d \. orig 
                        \/s_
                        (\d(?:_(?:1|2)){0,1})                           # lane and paired-end numbers
                        _sequence\.txt\.bz2

                        \s*
                        \z/xms ) { 
        $library_name  = $1;
        $flowcell_name = $2;
        $lane_end_nums = $3;

        if (! exists $labels{$library_name} ) { 
            die "Can't parse input: $input\n";
        }

        if ( $lane_end_nums =~ /\A (\d) _ (1|2) \z/xms ) { 
            $lane_num = $1;
            $end_num  = $2;
            if ( $end_num == 1) { 
                print 'ls -lt ', $prefix1A, $flowcell_name, $prefix1B, 's_', $lane_num, "_1_sequence.txt.bz2 ;\n";
                print 'ls -lt ', $prefix1A, $flowcell_name, $prefix1B, 's_', $lane_num, "_2_sequence.txt.bz2 ;\n";
                print 'bzcat ', $prefix1A, $flowcell_name, $prefix1B, 's_', $lane_num, "_1_sequence.txt.bz2 > $library_name.$flowcell_name.$lane_num.1.fq ;\n";
                print 'bzcat ', $prefix1A, $flowcell_name, $prefix1B, 's_', $lane_num, "_2_sequence.txt.bz2 > $library_name.$flowcell_name.$lane_num.2.fq ;\n";
                print "shuffleSequences_fastq.pl",
                      " $library_name.$flowcell_name.$lane_num.1.fq",
                      " $library_name.$flowcell_name.$lane_num.2.fq",
                      " $labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq ;",
                      "\n",
                      ;

                print "rm $library_name.$flowcell_name.$lane_num.1.fq $library_name.$flowcell_name.$lane_num.2.fq ;\n";
                print "gzip -9 $labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq ;\n";
                print "cp -ip $labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq.gz $directory1 ;\n";
            }
        }
        if ( $lane_end_nums =~ /\A (\d) \z/xms ) { 
            $lane_num = $1;
            $end_num  = q{};
            print 'ls -lt ', $prefix1A, $flowcell_name, $prefix1B, 's_', $lane_num, "_sequence.txt.bz2 ;\n";

            print 'bzcat ', $prefix1A, $flowcell_name, $prefix1B, 's_', $lane_num, "_sequence.txt.bz2 > $labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq ;\n";
            print "gzip -9 $labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq ;\n";
            print "cp -ip $labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq.gz $directory1 ;\n";
        }
        print "\n";
        print "# Note: tissue type $labels{$library_name}; library $library_name; flowcell $flowcell_name; lane $lane_num";
        if ($end_num) { 
            print "; end $end_num";
        }
        print "\n";
        print "# Visible at:  $prefix3", "$labels{$library_name}.$library_name.$flowcell_name.$lane_num.fq.gz \n"; 
    }
    else { 
        die "WARNING -- can't parse input line!: $input\n";
    }
}

