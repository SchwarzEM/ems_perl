#!/usr/bin/env perl

# build_job_06sep2011a.pl -- Erich Schwarz <emsch@caltech.edu>, 9/6/2011.
# Purpose: given a list of files, make a reliable script for handling their data.

use strict;
use warnings;

# Input: first_reads.06sep2011/2kb_Hco_02aug2011.paired.filt.noMOIdna.fq.gz.abundfilt.keep

my $stem   = q{};
my $suffix = q{};

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    if ( ( $input =~ /\A \S+ \z/xms ) and ( $input =~ /\A (\S+) \.paired\. (\S+) \z/xms ) ) { 
        $stem   = $1;
        $suffix = $2;

# Input: 2kb_Hco_02aug2011.paired.filt.noMOIdna.fq.gz.abundfilt.keep
# Stem:  2kb_Hco_02aug2011
#
#    filter_minimum_ACGT_fasta.pl -m 37 -i first_reads.06sep2011/2kb_Hco_02aug2011.paired.filt.noMOIdna.fq.gz.abundfilt.keep \
#        -o paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.pt1.txt ;
# 
#    filter_minimum_ACGT_fasta.pl -m 37 -i first_reads.06sep2011/2kb_Hco_02aug2011.unpaired.filt.noMOIdna.fq.gz.abundfilt.keep \
#        -o paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.pt2.txt ;

        my $text1  = '    filter_minimum_ACGT_fasta.pl -m 37 -i first_reads.06sep2011/';
        my $text3  = ' -o paired_reads.06sep2011/';
        my $text4a  = '.jumbled.pt1.txt';
        my $text4b  = '.jumbled.pt2.txt';

        my $int_file_1   = $stem . '.paired.' . $suffix;
        my $int_file_2   = $stem . '.unpaired.' . $suffix;

        my $int_file_3   = $stem . $text4a;
        my $int_file_4   = $stem . $text4b;

        my $int_file_5   = $stem . '.jumbled.06sep2011.filt1.fa';
        my $int_file_6   = $stem . '.jumbled.06sep2011.filt2.fa';

        my $out_paired   = $stem . '.paired.min37.noMOIdna2.06sep2011.fa';
        my $out_unpaired = $stem . '.unpaired.min37.noMOIdna2.06sep2011.fa';

        print "\n";

        print $text1, $int_file_1, $text3, $int_file_3, " ;\n", ;
        print "\n";

        print $text1, $int_file_2, $text3, $int_file_4, " ;\n", ;
        print "\n";

#    cat paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.pt1.txt paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.pt2.txt \
#        >  paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.06sep2011.filt1.fa ;

        my $text6  = '    cat paired_reads.06sep2011/';
        my $text7  = 'paired_reads.06sep2011/';

        print $text6, $int_file_3, q{ }, $text7, $int_file_4, q{ > }, $text7, $int_file_5, " ;\n", ;
        print "\n";

#    rm paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.pt1.txt paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.pt2.txt ;

        my $text8  = '    rm ';

        print $text8, $text7, $int_file_3, q{ }, $text7, $int_file_4, " ;\n", ;
        print "\n";

#     bowtie ../sheep.aug2011/Ovis_Bos_BTdb \
#        -p 12 -t --un paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.06sep2011.filt2.fa \
#        -f paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.06sep2011.filt1.fa /dev/null ;

        my $text9  = "    bowtie ../sheep.aug2011/Ovis_Bos_BTdb \\";
        my $text10 = '        -p 12 -t --un paired_reads.06sep2011/';
        my $text11 = '        -f paired_reads.06sep2011/';
        my $text12 = ' /dev/null ;';

        print $text9, "\n", ;
        print $text10, $int_file_6, " \\\n", ;
        print $text11, $int_file_5, $text12, " \n", ;
        print "\n";

#     paired_vs_unp_fastq.or.a.pl --fasta \
#        -i paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.06sep2011.filt2.fa \
#        -p paired_reads.06sep2011/2kb_Hco_02aug2011.paired.min37.noMOIdna2.06sep2011.fa \
#        -u unpaired_reads.06sep2011/2kb_Hco_02aug2011.unpaired.min37.noMOIdna2.06sep2011.fa ;

        my $text13 = "    paired_vs_unp_fastq.or.a.pl --fasta \\";
        my $text14 = '        -i paired_reads.06sep2011/';
        my $text15 = '        -p paired_reads.06sep2011/';
        my $text16 = '        -u unpaired_reads.06sep2011/';

        print $text13, " \n", ;
        print $text14, $int_file_6, " \\ \n", ;
        print $text15, $out_paired, " \\ \n", ;
        print $text16, $out_unpaired, " ; \n", ;
        print "\n";

#     rm paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.06sep2011.filt1.fa paired_reads.06sep2011/2kb_Hco_02aug2011.jumbled.06sep2011.filt2.fa ;

       print "    rm paired_reads.06sep2011/$int_file_5 paired_reads.06sep2011/$int_file_6 ;\n";
       print "\n";

    }
}

