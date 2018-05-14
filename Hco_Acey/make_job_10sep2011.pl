#!/usr/bin/env perl

# make_job_10sep2011.pl -- Erich Schwarz <emsch@caltech.edu>, 9/9/2011.
# Purpose: given a list of Titus-filtered files, generate an entire script for filtering, sorting, and assembling, with intermittent pings.

use strict;
use warnings;

my $in_file      = q{};
my $stem         = q{};
my $appendix     = q{};
my $int_file_1   = q{};
my $int_file_2   = q{};
my $int_file_3   = q{};
my $int_file_4   = q{};
my $out_paired   = q{};
my $out_unpaired = q{};

my %single_inputs  = ();
my %jumbled_inputs = ();
my %paired_inputs  = ();

my %final_unpaired = ();
my %final_paired   = ();

my @velvetg_sets = qw( shortPaired shortPaired2 shortPaired3 shortPaired4 shortPaired5 );

my %data2set = (   '42FC3AAXX.3'        => 'shortPaired',
                   '42YPHAAXX.6'        => 'shortPaired',
                   '61J8HAAXX.1'        => 'shortPaired',
                   '62636AAXX.1'        => 'shortPaired',

                   '625LPAAXX.4'        => 'shortPaired2',
                   '6263DAAXX.8'        => 'shortPaired2',
                   '62DJDAAXX.8'        => 'shortPaired2',
                   'B08VPABXX.3'        => 'shortPaired2',

                   '2kb_Hco_02aug2011'  => 'shortPaired3',

                   '5kb_Hco_02aug2011'  => 'shortPaired4',

                   '10kb_Hco_23jan2011' => 'shortPaired5',
);

while (my $input = <>) { 
    chomp $input;
    if ( ( $input !~ /\A \S+\.unpaired\.\S+ \z/xms ) and ( $input =~ /\A (\S+\.paired\.\S+) \z/xms ) ) { 
        $in_file = $1;
        $paired_inputs{$in_file} = 1;
    }
    elsif ( ( $input !~ /\A \S+\.unpaired\.\S+ \z/xms ) and ( $input =~ /\A (\S+\.jumbled\.\S+) \z/xms ) ) {
        $in_file = $1;
        $jumbled_inputs{$in_file} = 1;
    }
    elsif ( ( $input !~ /\A \S+\.unpaired\.\S+ \z/xms ) and ( $input =~ /\A (\S+\.single\.\S+) \z/xms ) ) {
        $in_file = $1;
        $single_inputs{$in_file} = 1;
    }
}

print '#!/bin/bash', "\n\n", ;

foreach my $single_infile (sort keys %single_inputs ) { 
    if ( $single_infile =~ /\A (\S+) \.single\.\S+ \z /xms ) { 
        $stem         = $1;
        $int_file_1   = $stem . '.single.10sep2011.filt1.fa';
        $out_unpaired = $stem . '.single.min35.noMOIdna2.10sep2011.fa';

        my $full_unpaired_name = '20th_filt_reads.09sep2011/' . $single_infile;
        if (! -e $full_unpaired_name ) {   
            die "Can't find: $full_unpaired_name\n";
        }

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i 20th_filt_reads.09sep2011/',
              $single_infile,
              ' -o unpaired_reads.10sep2011/',
              $int_file_1,
              " ;\n",
              ;

        print '    bowtie ../sheep.aug2011/Ovis_Bos_BTdb -p 12 -t --un unpaired_reads.10sep2011/',
              $out_unpaired,
              ' -f unpaired_reads.10sep2011/',
              $int_file_1,
              " /dev/null ;\n",
              ;

        print '    rm unpaired_reads.10sep2011/',
              $int_file_1,
              " ;\n",
              ;
        print "\n";

        my $final_file = 'unpaired_reads.10sep2011/' . $out_unpaired;
        if (exists $final_unpaired{$final_file}) { 
            die "Redundant final unpaired file: $final_file\n";
        }
        $final_unpaired{$final_file} = 1;
    }
    else { 
        die "Can't parse putative single-read infile: $single_infile\n";
    }
}

foreach my $jumbled_infile (sort keys %jumbled_inputs ) { 
    if ( $jumbled_infile =~ /\A (\S+) \.jumbled\.\S+ \z /xms ) { 
        $stem         = $1;
        $int_file_1   = $stem . '.jumbled.10sep2011.filt1.fa';
        $int_file_2   = $stem . '.jumbled.10sep2011.filt2.fa';
        $out_paired   = $stem . '.paired.min35.noMOIdna2.10sep2011.fa';
        $out_unpaired = $stem . '.unpaired.min35.noMOIdna2.10sep2011.fa';

        my $full_jumbled_name = '20th_filt_reads.09sep2011/' . $jumbled_infile;
        if (! -e $full_jumbled_name )  {
            die "Can't find: $full_jumbled_name\n";
        }

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i 20th_filt_reads.09sep2011/',
              $jumbled_infile,
              ' -o paired_reads.10sep2011/',
              $int_file_1,
              " ;\n",
              ;

        print '    bowtie ../sheep.aug2011/Ovis_Bos_BTdb -p 12 -t --un paired_reads.10sep2011/',
              $int_file_2,
              ' -f paired_reads.10sep2011/',
              $int_file_1,
              " /dev/null ;\n",
              ;

        print '    paired_vs_unp_fastq.or.a.pl --fasta',
              ' -i paired_reads.10sep2011/',
              $int_file_2,
              ' -p paired_reads.10sep2011/',
              $out_paired,
              ' -u unpaired_reads.10sep2011/',
              $out_unpaired,
              " ;\n",
              ;

        print '    rm paired_reads.10sep2011/',
              $int_file_1,
              ' paired_reads.10sep2011/',
              $int_file_2,
              " ;\n",
              ;
        print "\n";

        my $final_file = 'paired_reads.10sep2011/' . $out_paired;
        if (exists $final_paired{$final_file}) {
            die "Redundant final paired file: $final_file\n";
        }
        $final_paired{$final_file} = 1;
        
        $final_file = 'unpaired_reads.10sep2011/' . $out_unpaired;
        if (exists $final_unpaired{$final_file}) {
            die "Redundant final unpaired file: $final_file\n";
        }
        $final_unpaired{$final_file} = 1;
    }
    else { 
        die "Can't parse putative jumbled infile: $jumbled_infile\n";
    }
}

foreach my $paired_infile (sort keys %paired_inputs ) { 
    if ( $paired_infile =~ /\A (\S+) \.paired\. (\S+) \z /xms ) { 
        $stem         = $1;
        $appendix     = $2;
        my $unpaired  = $stem . '.unpaired.' . $appendix;
        $int_file_1   = $stem . '.jumbled.10sep2011.filt1.fa';
        $int_file_2   = $stem . '.jumbled.10sep2011.filt2.fa';
        $int_file_3   = $stem . '.jumbled.10sep2011.filt3.fa';
        $int_file_4   = $stem . '.jumbled.10sep2011.filt4.fa';
        $out_paired   = $stem . '.paired.min35.noMOIdna2.10sep2011.fa';
        $out_unpaired = $stem . '.unpaired.min35.noMOIdna2.10sep2011.fa';

        my $full_paired_name   = '20th_filt_reads.09sep2011/' . $paired_infile;
        my $full_unpaired_name = '20th_filt_reads.09sep2011/' . $unpaired;
        if (! -e $full_paired_name )  {   
            die "Can't find: $full_paired_name\n";
        }
        if (! -e $full_unpaired_name ) {
            die "Can't find: $full_unpaired_name\n";
        }

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i 20th_filt_reads.09sep2011/',
              $paired_infile,
              ' -o paired_reads.10sep2011/',
              $int_file_1,
              " ;\n",
              ;

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i 20th_filt_reads.09sep2011/',
              $unpaired,
              ' -o paired_reads.10sep2011/',
              $int_file_2,
              " ;\n",
              ;

        print '    cat paired_reads.10sep2011/',
              $int_file_1,
              ' paired_reads.10sep2011/',
              $int_file_2,
              ' > paired_reads.10sep2011/',
              $int_file_3,
              " ; \n",
              ;

        print '    bowtie ../sheep.aug2011/Ovis_Bos_BTdb -p 12 -t --un paired_reads.10sep2011/',
              $int_file_4,
              ' -f paired_reads.10sep2011/',
              $int_file_3,
              " /dev/null ;\n",
              ;

        print '    paired_vs_unp_fastq.or.a.pl --fasta',
              ' -i paired_reads.10sep2011/',
              $int_file_4,
              ' -p paired_reads.10sep2011/',
              $out_paired,
              ' -u unpaired_reads.10sep2011/',
              $out_unpaired,
              " ;\n",
              ;

        print '    rm paired_reads.10sep2011/',
              $int_file_1,
              ' paired_reads.10sep2011/',
              $int_file_2,
              ' paired_reads.10sep2011/',
              $int_file_3,
              ' paired_reads.10sep2011/',
              $int_file_4,
              " ;\n",
              ;
        print "\n";

        my $final_file = 'paired_reads.10sep2011/' . $out_paired;
        if (exists $final_paired{$final_file}) {
            die "Redundant final paired file: $final_file\n";
        }
        $final_paired{$final_file} = 1;
        
        $final_file = 'unpaired_reads.10sep2011/' . $out_unpaired;
        if (exists $final_unpaired{$final_file}) {
            die "Redundant final unpaired file: $final_file\n";
        }
        $final_unpaired{$final_file} = 1;
    }
    else { 
        die "Can't parse putative jumbled paired-read infile: $paired_infile\n";
    }
}

print "    program_done_e-ping.pl -p filtering_sorting_reads_10sep2011 ;\n";
print "\n";

my @unpaired_finals = sort keys %final_unpaired;
my @paired_finals   = sort keys %final_paired;

print "    velveth Hco_VelGen_10sep2011_allBGI_standard_min35_k35 35 \\ \n";
print "        -fasta \\ \n";
print "        -short \\ \n";

foreach my $unp1 (@unpaired_finals) { 
    print "        $unp1 \\ \n";
}

print "                  \\ \n";

foreach my $vg_set (@velvetg_sets) { 
    print '    -', 
          "$vg_set  \\ \n", 
          ;
    foreach my $pai1 (@paired_finals) {
        if ( $pai1 =~ /\A paired_reads.10sep2011\/ (\S+) \.paired\. \S+ \z/xms ) { 
            $stem = $1;
            if (! exists $data2set{$stem} ) { 
                die "Can't assign $pai1 to velvetg group on basis of $stem\n";
            }
            if ( $data2set{$stem} eq $vg_set ) { 
                print "        $pai1 \\ \n";
            }
        }
        else { 
            die "Can't parse putative paired file: $pai1\n";
        }
    }
}

print "    ; \n";
print "\n";
print "    program_done_e-ping.pl -p velveth_10sep2011_allBGI_min35_k35 ; \n";
print "\n";
print "    velvetg Hco_VelGen_10sep2011_allBGI_standard_min35_k35 \\ \n";
print "        -shortMatePaired3 yes -shortMatePaired4 yes -shortMatePaired5 yes \\ \n";
print "        -cov_cutoff 4 -exp_cov 100 -min_contig_lgth 200 \\ \n";
print "        -ins_length  300 -ins_length_sd  50 \\ \n";
print "        -ins_length2 500 -ins_length2_sd 200 \\ \n";
print "        -ins_length3  2000 \\ \n";
print "        -ins_length4  5000 \\ \n";
print "        -ins_length5 10000 ; \n";
print "\n";
print "    program_done_e-ping.pl -p velvet_assembly_10sep2011_allBGI_min35_k35 ; \n";
print "\n";
print "    count_fasta_residues.pl Hco_VelGen_10sep2011_allBGI_standard_min35_k35/contigs.fa > full.count.10sep2011.txt ; \n";
print "    fasta_size_subset.pl -s --max 315M -i Hco_VelGen_10sep2011_allBGI_standard_min35_k35/contigs.fa | count_fasta_residues.pl > top.count.10sep2011.txt ; \n";
print "\n";

