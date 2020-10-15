#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};

$infile = $ARGV[0] if $ARGV[0];
open my $INFILE, '<', $infile;

my $out_porechop = $infile . '.porechop';
$out_porechop = safename($out_porechop);
open my $PORECHOP, '>', $out_porechop;

my $out_racon = $infile . '.racon';
$out_racon = safename($out_racon);
open my $RACON, '>', $out_racon;

my $out_medaka = $infile . '.medaka';
$out_medaka = safename($out_medaka);
open my $MEDAKA, '>', $out_medaka;

my $out_pilon = $infile . '.pilon';   # pilon|pil_
$out_pilon = safename($out_pilon);
open my $PILON, '>', $out_pilon;

my $out_ngm = $infile . '.ngm';
$out_ngm = safename($out_ngm);
open my $NGM, '>', $out_ngm;

my $out_ntEdit = $infile . '.ntEdit';
$out_ntEdit = safename($out_ntEdit);
open my $NTEDIT, '>', $out_ntEdit;

my $out_polca = $infile . '.polca';
$out_polca = safename($out_polca);
open my $POLCA, '>', $out_polca;

my $out_sourmash = $infile . '.sourmash';
$out_sourmash = safename($out_sourmash);
open my $SOURMASH, '>', $out_sourmash;

my $out_purge_dups = $infile . '.purge_dups';
$out_purge_dups = safename($out_purge_dups);
open my $PURGEDUPS, '>', $out_purge_dups;

my $out_minimap2 = $infile . '.minimap2';
$out_minimap2 = safename($out_minimap2);
open my $MINIMAP2, '>', $out_minimap2;

my $out_fastp = $infile . '.fastp';
$out_fastp = safename($out_fastp);
open my $FASTP, '>', $out_fastp;

my $out_hisat2 = $infile . '.hisat2';
$out_hisat2 = safename($out_hisat2);
open my $HISAT2, '>', $out_hisat2;

my $out_bowtie2 = $infile . '.bowtie2';  # bowtie2|bt2
$out_bowtie2 = safename($out_bowtie2);
open my $BOWTIE2, '>', $out_bowtie2;

my $out_bwa = $infile . '.bwa';
$out_bwa = safename($out_bwa);
open my    $BWA, '>', $out_bwa;

my $out_quickmerge = $infile . '.quickmerge';
$out_quickmerge = safename($out_quickmerge);
open my $QUICKMERGE, '>', $out_quickmerge;

my $out_LRScaf = $infile . '.LRScaf';
$out_LRScaf = safename($out_LRScaf);
open my $LRSCAF, '>', $out_LRScaf;

my $out_rmask = $infile . '.rmask';
$out_rmask = safename($out_rmask);
open my $RMASK, '>', $out_rmask;

my $out_busco = $infile . '.busco';
$out_busco = safename($out_busco);
open my $BUSCO, '>', $out_busco;

my $out_blast = $infile . '.blast';
$out_blast = safename($out_blast);
open my $BLAST, '>', $out_blast;

my $out_orthofinder = $infile . '.orthofinder';
$out_orthofinder = safename($out_orthofinder);
open my $OFINDER, '>', $out_orthofinder;

my $out_pfam = $infile . '.Pfam';  # Pfam|hmms
$out_pfam = safename($out_pfam);
open my $PFAM, '>', $out_pfam;

my $out_interpro = $infile . '.interpro';  # InterPro|interpro|PantherDB
$out_interpro = safename($out_interpro);
open my $INTERPRO, '>', $out_interpro;

my $out_phobius = $infile . '.phobius';
$out_phobius = safename($out_phobius);
open my $PHOBIUS, '>', $out_phobius;

my $out_khmer = $infile . '.khmer';
$out_khmer = safename($out_khmer);
open my $KHMER, '>', $out_khmer;

my $out_braker2 = $infile . '.braker2';   # braker2|augustus
$out_braker2 = safename($out_braker2);
open my $BRAKER2, '>', $out_braker2;

my $out_salmon = $infile . '.salmon';
$out_salmon = safename($out_salmon);
open my $SALMON, '>', $out_salmon;

my $out_wrangle = $infile . '.wrangle';
$out_wrangle = safename($out_wrangle);
open my $WRANGLE, '>', $out_wrangle;

my $out_canu = $infile . '.canu';  # Canu|canu
$out_canu = safename($out_canu);
open my $CANU, '>', $out_canu;

my $out_soapdenovo = $infile . '.soapdenovo';
$out_soapdenovo = safename($out_soapdenovo);
open my $SOAPDENOVO, '>', $out_soapdenovo;

my $out_flye = $infile . '.flye';   # _Flye_|_flye_
$out_flye = safename($out_flye);
open my $FLYE, '>', $out_flye;

my $out_raven = $infile . '.raven';
$out_raven = safename($out_raven);
open my $RAVEN, '>', $out_raven;

my $out_wtdbg2 = $infile . '.wtdbg2';
$out_wtdbg2 = safename($out_wtdbg2);
open my $WTDBG2, '>', $out_wtdbg2;

my $out_unclassified = $infile . '.unclassified';
$out_unclassified = safename($out_unclassified);
open my  $UNCLASS, '>', $out_unclassified;

while (my $input = <$INFILE>) {
    chomp $input;
    $input =~ s/\s+\z//;

    # Sample input line:
    # 10820792    job_raxei_Eur_2020.08.10.phobius_2020.09.03.01.sh    2020-09-10 18:31:18    \
    # 2020-09-10 19:58:54    RM-SHARED    1    87.6    1.46

    if ( $input =~ /\A [^\t]* \t (\S+) \t .* \t (\S+) \z/xms ) {
        my $job_name = $1;
        my $cpus     = $2;

        if ( $job_name =~ /porechop/xms ) {
             print $PORECHOP "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /racon/xms ) {
             print $RACON "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /medaka/xms ) {
             print $MEDAKA "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /pilon|pil_/xms  ) {
             print $PILON "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /ngm/xms ) {
             print $NGM "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /ntEdit|ntHit/xms ) {
             print $NTEDIT "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /polca/xms ) {
             print $POLCA "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /sourmash/xms ) {
             print $SOURMASH "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /purge_dups/xms ) {
             print $PURGEDUPS "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /minimap2/xms ) {
             print $MINIMAP2 "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /fastp/xms ) {
             print $FASTP "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /hisat2/xms ) {
             print $HISAT2 "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /bowtie2|bt2/xms ) {
             print $BOWTIE2 "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /bwa/xms ) {
             print $BWA "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /quickmerge/xms ) {
             print $QUICKMERGE "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /LRScaf|rnascaf/xms ) {
             print $LRSCAF "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /rmask|reps|RepeatModeler/xms ) {
             print $RMASK "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /busco/xms ) {
             print $BUSCO "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /blast/xms  ) {
             print $BLAST "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /orthofinder/xms ) {
             print $OFINDER "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /Pfam|hmms/xms ) {
             print $PFAM "$job_name\t$cpus\n";
        }
	elsif ( $job_name =~ /InterPro|interpro|PantherDB/xms ) {
             print $INTERPRO "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /phobius/xms ) {
             print $PHOBIUS "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /khmer/xms ) {
             print $KHMER "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /braker2|augustus|bam2hints/xms ) {
             print $BRAKER2 "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /salmon|gentrome/xms ) {
             print $SALMON "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /wrangl|gzip|pool|split20|unpack|wget|unzip/xms ) {
             print $WRANGLE "$job_name\t$cpus\n";
        }
	elsif ( $job_name =~ /Canu|canu/xms ) {
             print $CANU "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /soapdenovo/xms ) {
             print $SOAPDENOVO "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /_Flye_|_flye_/xms ) {
             print $FLYE "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /_raven_/xms ) {
             print $RAVEN "$job_name\t$cpus\n";
        }
        elsif ( $job_name =~ /_wtdbg2_/xms ) {
             print $WTDBG2 "$job_name\t$cpus\n";
        }
        else {
             print $UNCLASS "$job_name\t$cpus\n";
        }
    }
    elsif ( $input !~ /\AJob ID/ ) {
        die "Cannot parse input line from $infile: $input\n";
    }
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

