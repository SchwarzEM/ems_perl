#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw(uniq);
use Getopt::Long;

my $geneids  = q{};
my @phenos   = ();
my @nophenos = ();
my $help;

my $header = "Gene\tWBPheno\tWBNonPheno";

my $data_ref;

GetOptions ( 'geneids=s'     => \$geneids,
             'phenos=s{,}'   => \@phenos,
             'nophenos=s{,}' => \@nophenos,
             'help'          => \$help,   );

if ( $help or (! $geneids) or (! @phenos) or (! @nophenos) ) {
    die "Format: gene2pheno_table_03feb2015.pl\n",
        "    --geneids|-g    <single WormBase gene ID file from (e.g.) WS245>\n",
        "    --phenos|-p     <one or more files with WormBase gene phenotypes from RNAi or alleles,>\n",
        "    --nophenos|-n   <one or more files with WormBase gene no-phenotypes from RNAi or alleles>\n",
        "                    [note that all three tables should ideally be from the same WormBase release]\n",
        "    --help|-h       [print this message]\n",
        ;
}

foreach my $pheno_file (@phenos) {
    open my $PHENO, '<', $pheno_file;
    while (my $input = <$PHENO>) {
        chomp $input;

        # Sample inputs:
        # "WBVar00000005"  "WBGene00001134"  "WBPhenotype:0000640"  "egg laying variant"    "oviposition abnormal"
        # "WBVar00000005"  "WBGene00001134"  "WBPhenotype:0000646"  "sluggish"      "lethargic"
        # "WBVar00000005"  "WBGene00001134"  "WBPhenotype:0000646"  "sluggish"      "slow"
        # "WBVar00000005"  "WBGene00001134"  "WBPhenotype:0000646"  "sluggish"      "Slu"
        # "WBVar00000005"  "WBGene00001134"  "WBPhenotype:0000886"  "Variant"       "Abnormal"
        # "WBVar00000001"  "WBGene00003883"  "WBPhenotype:0000843"  "male mating efficiency reduced"     * no synonym column *

        if ( $input =~ /\A \"WB[A-Za-z]{3,4}\d+\" \t \"(WBGene\d+)\" \t \"(WBPhenotype:\d+)\" \t ([^\t]+) \t ([^\t]*) \z/xms ) {
            my $wbgene     = $1;
            my $wbpheno    = $2;
            my $pheno_desc = $3;
            my $pheno_syn  = $4;

            # Strip off leading and trailing double-quotes:
            $pheno_desc =~ s/\A\"//;
            $pheno_desc =~ s/\"\z//;
            $pheno_syn  =~ s/\A\"//; 
            $pheno_syn  =~ s/\"\z//;

            $data_ref->{'wbgene'}->{$wbgene}->{'wbpheno'}->{$wbpheno} = 1;

            if ( ( exists $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} ) and ( $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} ne $pheno_desc ) ) {
                die "For $wbpheno, inconsistent descriptions: $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} versus $pheno_desc\n";
            }

            $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} = $pheno_desc;
            if ($pheno_syn) {
                $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_syn'}->{$pheno_syn} = 1;
            }
        }
        else { 
            die "From phenotype file $pheno_file, can't parse: $input\n";
        }
    }
    close $PHENO;
}

foreach my $no_pheno_file (@nophenos) {
    open my $NOPHENO, '<', $no_pheno_file;
    while (my $input = <$NOPHENO>) {
        chomp $input;
            
        # Sample inputs:
        # "WBVar00000005" "WBGene00001134"        "WBPhenotype:0000640"   "egg laying variant"    "oviposition abnormal"
        # "WBVar00000005" "WBGene00001134"        "WBPhenotype:0000646"   "sluggish"      "lethargic"
        # "WBVar00000005" "WBGene00001134"        "WBPhenotype:0000646"   "sluggish"      "slow"
        # "WBVar00000005" "WBGene00001134"        "WBPhenotype:0000646"   "sluggish"      "Slu"
        # "WBVar00000005" "WBGene00001134"        "WBPhenotype:0000886"   "Variant"       "Abnormal"

        if ( $input =~ /\A \"WB[A-Za-z]{3,4}\d+\" \t \"(WBGene\d+)\" \t \"(WBPhenotype:\d+)\" \t ([^\t]+) \t ([^\t]*) \z/xms ) {
            my $wbgene     = $1;
            my $wbpheno    = $2;
            my $pheno_desc = $3;
            my $pheno_syn  = $4;

            # Strip off leading and trailing double-quotes:
            $pheno_desc =~ s/\A\"//;
            $pheno_desc =~ s/\"\z//; 
            $pheno_syn  =~ s/\A\"//;
            $pheno_syn  =~ s/\"\z//;

            $data_ref->{'wbgene'}->{$wbgene}->{'no_wbpheno'}->{$wbpheno} = 1;

            if ( ( exists $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} ) and ( $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} ne $pheno_desc ) ) {
                die "For $wbpheno, inconsistent descriptions: $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} versus $pheno_desc\n";
            }
        
            $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'} = $pheno_desc;
            if ($pheno_syn) {
                $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_syn'}->{$pheno_syn} = 1;
            }
        }
        else {
            die "From no-phenotype file $no_pheno_file, can't parse: $input\n";
        }
    }
    close $NOPHENO;
}

# We do not want to list IDs for genes for which we have no data, so the ID table gets parsed last.

open my $GENEIDS, '<', $geneids;
while (my $input = <$GENEIDS>) {
    chomp $input;
        
    # Sample input:
    # WBGene00000001  aap-1   Y110A7A.10
    # WBGene00000002  aat-1   F27C8.1
    # [...]
    # WBGene00001590  C45G3.3 C45G3.3
    # WBGene00001608  R07B1.8 R07B1.8
     
    my $geneid = q{};

    if ( $input =~ /\A (WBGene\d+) \t (\S+) \t (\S+) \z/xms ) {
        my $wbgene = $1;
        my $cgc    = $2;
        my $cds    = $3;
        if ( $cds ne $cgc ) {
            $geneid = $wbgene . q{|} . $cds . q{|} . $cgc;
        }
        else {
            $geneid = $wbgene . q{|} . $cds;
        }
        if ( ( exists $data_ref->{'wbgene'}->{$wbgene}->{'geneid'} ) and ( $data_ref->{'wbgene'}->{$wbgene}->{'geneid'} ne $geneid ) ) {
            die "For WBGene ID $wbgene, two inconsistent full names: $data_ref->{'wbgene'}->{$wbgene}->{'geneid'} and $geneid\n";
        }
        if ( exists $data_ref->{'wbgene'}->{$wbgene} ) {
            $data_ref->{'wbgene'}->{$wbgene}->{'geneid'} = $geneid;
        }
    }
    else {
        die "From gene ID table $geneids, can't parse: $input\n";
    }
}       
close $GENEIDS;

my @wbgenes = sort keys %{ $data_ref->{'wbgene'} };

foreach my $wbgene (@wbgenes) {
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'geneid'} ) { 
        my $geneid = $data_ref->{'wbgene'}->{$wbgene}->{'geneid'};

        my @pheno_descs         = ();
        my @wbphenos            = ();
        my $full_pheno_desc_txt = q{};

        my @no_pheno_descs         = ();
        my @no_wbphenos            = ();
        my $full_no_pheno_desc_txt = q{};

        if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'wbpheno'} ) {
            @wbphenos = sort keys %{ $data_ref->{'wbgene'}->{$wbgene}->{'wbpheno'} };

            foreach my $wbpheno (@wbphenos) {
                my $pheno_desc    = $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_desc'};
                my @pheno_syns    = ();
                my $pheno_syn_txt = q{};

                if ( exists $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_syn'} ) {
                    @pheno_syns = sort keys %{ $data_ref->{'wbpheno'}->{$wbpheno}->{'pheno_syn'} };
                    $pheno_syn_txt  = join ', ', @pheno_syns;
                }

                my $pheno_desc_txt = $pheno_desc;
                if ($pheno_syn_txt) {
                    $pheno_desc_txt = $pheno_desc . " [$pheno_syn_txt]";
                }
                push @pheno_descs, $pheno_desc_txt;
            }

            @pheno_descs = sort @pheno_descs;
            @pheno_descs = uniq @pheno_descs;
            $full_pheno_desc_txt = join '; ', @pheno_descs;
        }

        if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'no_wbpheno'} ) {
            @no_wbphenos = sort keys %{ $data_ref->{'wbgene'}->{$wbgene}->{'no_wbpheno'} };

            foreach my $no_wbpheno (@no_wbphenos) {    
                my $pheno_desc    = $data_ref->{'wbpheno'}->{$no_wbpheno}->{'pheno_desc'};
                my @pheno_syns    = ();
                my $pheno_syn_txt = q{};
            
                if ( exists $data_ref->{'wbpheno'}->{$no_wbpheno}->{'pheno_syn'} ) {
                    @pheno_syns = sort keys %{ $data_ref->{'wbpheno'}->{$no_wbpheno}->{'pheno_syn'} };
                    $pheno_syn_txt  = join ', ', @pheno_syns;
                }

                my $pheno_desc_txt = $pheno_desc;
                if ($pheno_syn_txt) {
                    $pheno_desc_txt = $pheno_desc . " [$pheno_syn_txt]"; 
                } 
                push @no_pheno_descs, $pheno_desc_txt;
            }
            @no_pheno_descs = sort @no_pheno_descs;
            @no_pheno_descs = uniq @no_pheno_descs;
            $full_no_pheno_desc_txt = join '; ', @no_pheno_descs;
        }

        if ( @wbphenos or @no_wbphenos ) {
            if ($full_pheno_desc_txt) {
                $full_pheno_desc_txt = "WBPheno: $full_pheno_desc_txt";
            }
            if ($full_no_pheno_desc_txt) {
                $full_no_pheno_desc_txt = "WBNonPheno: $full_no_pheno_desc_txt";
            }

            # We only one this once, at the start of the file:
            print "$header\n" if $header;
            $header = q{};

            # Data line output:
            print "$geneid\t$full_pheno_desc_txt\t$full_no_pheno_desc_txt\n";
        }
    }
}

