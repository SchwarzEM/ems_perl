#!/usr/bin/env perl

# make_job_11sep2011.pl -- Erich Schwarz <emsch@caltech.edu>, 9/11/2011.
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

my $older_date       = '10sep2011';
my $newer_date       = '11sep2011';

my $in_reads_dir     = 'khmer_filt_reads.' . $older_date;
my $out_unpaired_dir = 'unpaired_reads.'   . $newer_date;
my $out_paired_dir   = 'paired_reads.'     . $newer_date;

my %single_inputs  = ();
my %jumbled_inputs = ();
my %paired_inputs  = ();

my %final_unpaired = ();
my %final_paired   = ();

my @velvetg_sets = qw( shortPaired shortPaired2 shortPaired3 shortPaired4 );

my %data2set = (   '11080_42FC3AAXX_c152_l8_r12'  => 'shortPaired',    # 11080      395 +/- 35 nt   42FC3AAXX, lane 8

                   '11079_42FC3AAXX_c152_l7_r12'  => 'shortPaired2',    # 11079      440 +/- 45 nt   42FC3AAXX, lane 7

                   '10513_30DY0AAXX_c151_l1_r12'  => 'shortPaired3',    # 10513      460 +/- 80 nt   10513_30DY0AAXX_c151_l1_r12,
                   '10513_ilmn200901_c202_l1_r12' => 'shortPaired3',    #                                10513_ilmn200901_c202_l1_r12,
                   'B08VPABXX.1'                  => 'shortPaired3',    #                                and B08VPABXX, lane 1

                   '616KLAAXX.1_jump'             => 'shortPaired4',    # needs empirical characterization, but could be ~3000 nt
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
        $int_file_1   = $stem . '.single.' . $newer_date . '.filt1.fa';
        $out_unpaired = $stem . '.single.min35.noMOIdna2.' . $newer_date . '.fa';

        my $full_unpaired_name = $in_reads_dir . q{/} . $single_infile;
        if (! -e $full_unpaired_name ) {   
            warn "Can't find: $full_unpaired_name\n";
        }

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i ', $in_reads_dir, q{/},
              $single_infile,
              ' -o ', $out_unpaired_dir, q{/}, 
              $int_file_1,
              " ;\n",
              ;

        print '    bowtie ../sheep.aug2011/Ovis_Bos_BTdb -p 12 -t --un ', $out_unpaired_dir, q{/}, 
              $out_unpaired,
              ' -f ', $out_unpaired_dir, q{/},
              $int_file_1,
              " /dev/null ;\n",
              ;

        print '    rm ', $out_unpaired_dir, q{/}, 
              $int_file_1,
              " ;\n",
              ;
        print "\n";

        my $final_file =  $out_unpaired_dir . q{/} . $out_unpaired;
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
        $int_file_1   = $stem . '.jumbled.' . $newer_date . '.filt1.fa';
        $int_file_2   = $stem . '.jumbled.' . $newer_date . '.filt2.fa';
        $out_paired   = $stem . '.paired.min35.noMOIdna2.' . $newer_date . '.fa';
        $out_unpaired = $stem . '.unpaired.min35.noMOIdna2.' . $newer_date . '.fa';

        my $full_jumbled_name =  $in_reads_dir . q{/} . $jumbled_infile;
        if (! -e $full_jumbled_name )  {
            warn "Can't find: $full_jumbled_name\n";
        }

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i ', $in_reads_dir, q{/},
              $jumbled_infile,
              ' -o ', $out_paired_dir, q{/}, 
              $int_file_1,
              " ;\n",
              ;

        print '    bowtie ../sheep.aug2011/Ovis_Bos_BTdb -p 12 -t --un ', $out_paired_dir, q{/}, 
              $int_file_2,
              ' -f ', $out_paired_dir, q{/}, 
              $int_file_1,
              " /dev/null ;\n",
              ;

        print '    paired_vs_unp_fastq.or.a.pl --fasta',
              ' -i ', $out_paired_dir, q{/}, 
              $int_file_2,
              ' -p ', $out_paired_dir, q{/}, 
              $out_paired,
              ' -u ', $out_unpaired_dir, q{/}, 
              $out_unpaired,
              " ;\n",
              ;

        print '    rm ', $out_paired_dir, q{/}, 
              $int_file_1,
              ' ', $out_paired_dir, q{/}, 
              $int_file_2,
              " ;\n",
              ;
        print "\n";

        my $final_file = $out_paired_dir . q{/} . $out_paired;
        if (exists $final_paired{$final_file}) {
            die "Redundant final paired file: $final_file\n";
        }
        $final_paired{$final_file} = 1;
        
        $final_file = $out_unpaired_dir . q{/} . $out_unpaired;
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
        $int_file_1   = $stem . '.jumbled.' . $newer_date . '.filt1.fa';
        $int_file_2   = $stem . '.jumbled.' . $newer_date . '.filt2.fa';
        $int_file_3   = $stem . '.jumbled.' . $newer_date . '.filt3.fa';
        $int_file_4   = $stem . '.jumbled.' . $newer_date . '.filt4.fa';
        $out_paired   = $stem . '.paired.min35.noMOIdna2.' . $newer_date . '.fa';
        $out_unpaired = $stem . '.unpaired.min35.noMOIdna2.' . $newer_date . '.fa';

        my $full_paired_name   = $in_reads_dir . q{/} . $paired_infile;
        my $full_unpaired_name = $in_reads_dir . q{/} . $unpaired;
        if (! -e $full_paired_name )  {   
            warn "Can't find: $full_paired_name\n";
        }
        if (! -e $full_unpaired_name ) {
            warn "Can't find: $full_unpaired_name\n";
        }

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i ', $in_reads_dir, q{/},
              $paired_infile,
              ' -o ', $out_paired_dir, q{/}, 
              $int_file_1,
              " ;\n",
              ;

        print '    filter_minimum_ACGT_fasta.pl -m 35 -i ', $in_reads_dir, q{/},
              $unpaired,
              ' -o ', $out_paired_dir, q{/}, 
              $int_file_2,
              " ;\n",
              ;

        print '    cat ', $out_paired_dir, q{/}, 
              $int_file_1,
              ' ', $out_paired_dir, q{/}, 
              $int_file_2,
              ' > ', $out_paired_dir, q{/}, 
              $int_file_3,
              " ; \n",
              ;

        print '    bowtie ../sheep.aug2011/Ovis_Bos_BTdb -p 12 -t --un ', $out_paired_dir, q{/}, 
              $int_file_4,
              ' -f ', $out_paired_dir, q{/}, 
              $int_file_3,
              " /dev/null ;\n",
              ;

        print '    paired_vs_unp_fastq.or.a.pl --fasta',
              ' -i ', $out_paired_dir, q{/}, 
              $int_file_4,
              ' -p ', $out_paired_dir, q{/}, 
              $out_paired,
              ' -u ', $out_unpaired_dir, q{/}, 
              $out_unpaired,
              " ;\n",
              ;

        print '    rm ', $out_paired_dir, q{/},
              $int_file_1,
              ' ', $out_paired_dir, q{/}, 
              $int_file_2,
              ' ', $out_paired_dir, q{/}, 
              $int_file_3,
              ' ', $out_paired_dir, q{/}, 
              $int_file_4,
              " ;\n",
              ;
        print "\n";

        my $final_file = $out_paired_dir . q{/} . $out_paired;
        if (exists $final_paired{$final_file}) {
            die "Redundant final paired file: $final_file\n";
        }
        $final_paired{$final_file} = 1;
        
        $final_file = $out_unpaired_dir . q{/} . $out_unpaired;
        if (exists $final_unpaired{$final_file}) {
            die "Redundant final unpaired file: $final_file\n";
        }
        $final_unpaired{$final_file} = 1;
    }
    else { 
        die "Can't parse putative jumbled paired-read infile: $paired_infile\n";
    }
}

print "    program_done_e-ping.pl -p filtering_sorting_reads_", "$newer_date ;\n";
print "\n";

my @unpaired_finals = sort keys %final_unpaired;
my @paired_finals   = sort keys %final_paired;

print "    velveth Hco_VelGen_", $newer_date, "_allBGI_standard_min35_k35 35 \\ \n";
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
        if ( $pai1 =~ /\A $out_paired_dir \/ (\S+) \.paired\. \S+ \z/xms ) { 
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

print "    program_done_e-ping.pl -p velveth_", $newer_date, "_allBGI_min35_k35 ; \n";
print "\n";

print "    velvetg Hco_VelGen_", $newer_date, "_allBGI_standard_min35_k35 \\ \n";
print "        -shortMatePaired4 yes \\ \n";
print "        -cov_cutoff 4 -exp_cov 100 -min_contig_lgth 200 \\ \n";
print "        -ins_length  395 -ins_length_sd  35 \\ \n";
print "        -ins_length2 440 -ins_length2_sd 45 \\ \n";
print "        -ins_length3 460 -ins_length3_sd 80 \\ \n";
print "        -ins_length4 3000 ; \n";
print "\n";

print "    program_done_e-ping.pl -p velvet_assembly_", $newer_date, "_allBGI_min35_k35 ; \n";
print "\n";

print "    count_fasta_residues.pl Hco_VelGen_", 
      $newer_date, 
      "_allBGI_standard_min35_k35/contigs.fa > full.count.", 
      $newer_date, 
      ".txt ; \n",
      ;

print "    fasta_size_subset.pl -s --max 315M -i Hco_VelGen_", 
      $newer_date, 
      "_allBGI_standard_min35_k35/contigs.fa | count_fasta_residues.pl > top.count.", 
      $newer_date, ".txt ; \n",
      ;
print "\n";

