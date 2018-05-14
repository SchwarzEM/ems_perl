#!/usr/bin/env perl

use strict;
use warnings;

my $full_tfs = $ARGV[0];

my $tfs_whf  = $ARGV[1];
my $tfs_wh   = $ARGV[2];
my $tfs_wf   = $ARGV[3];

my $tfs_whf2 = $ARGV[4];
my $tfs_wh2  = $ARGV[5];

my $tfs_reg      = $ARGV[6];
my $tfs_flk      = $ARGV[7];
my $tfs_cons_flk = $ARGV[8];

my @known_cons    = ($tfs_whf, $tfs_wh, $tfs_wf);
my @possible_cons = ($tfs_whf2, $tfs_wh2);

my $data_ref;

open my $FULL_TFS, '<', $full_tfs or die "Can't open full TF annotation file, $full_tfs: $!";
while (my $input = <$FULL_TFS> ) { 
    chomp $input;
    # Sample input line:
    # WBGene00019424  K06A1.1 aptf-1  AP-2    ENSP00000263543 TFAP2E  "1, 2"
    if ( $input =~ /\A (WBGene\d+) \t [^\t]* \t [^\t]* \t ([^\t]+) /xms ) { 
        my $wbgene  = $1;
        my $tf_type = $2;
        $data_ref->{'wbgene'}->{$wbgene}->{'tf_type'} = $tf_type;
        $data_ref->{'key_wbgene'}->{$wbgene} = 1;
    }
    elsif ( $input =~ /\A (WBGene\d+) /xms ) { 
        die "From full TF annotation file, $full_tfs, can't parse input line: $input\n";
    }
}
close $FULL_TFS or die "Can't close filehandle to full TF annotation file, $full_tfs: $!";

foreach my $known_tf_cons (@known_cons) { 
    open my $KNOWN_TF_CONS, '<', $known_tf_cons or die "Can't open TF known conservation file, $known_tf_cons: $!";
    while (my $input = <$KNOWN_TF_CONS> ) {
        chomp $input;
        if ( $input =~ /\A (WBGene\d+) \t (?: [^\t]* \t){6} ([^\t]+) /xms ) {
            my $wbgene  = $1;
            my $treefam = $2;
            $data_ref->{'wbgene'}->{$wbgene}->{'known_treefam'} = $treefam;
        }
        elsif ( $input =~ /\A (WBGene\d+) /xms ) { 
            die "From TF known conservation file, $known_tf_cons, can't parse input line: $input\n";
        }
    }
    close $KNOWN_TF_CONS or die "Can't close filehandle to TF known conservation file, $known_tf_cons: $!";
}

foreach my $possible_tf_cons (@possible_cons) {
    open my $POSSIBLE_TF_CONS, '<', $possible_tf_cons or die "Can't open TF possible conservation file, $possible_tf_cons: $!";
    while (my $input = <$POSSIBLE_TF_CONS> ) {
        chomp $input;
        if ( $input =~ /\A (WBGene\d+) \t (?: [^\t]* \t){6} ([^\t]+) /xms ) {
            my $wbgene  = $1;
            my $treefam = $2;
            $data_ref->{'wbgene'}->{$wbgene}->{'possible_treefam'} = $treefam;
        }
        elsif ( $input =~ /\A (WBGene\d+) /xms ) {
            die "From TF possible conservation file, $possible_tf_cons, can't parse input line: $input\n";
        }
    }
    close $POSSIBLE_TF_CONS or die "Can't close filehandle to TF possible conservation file, $possible_tf_cons: $!";
}

open my $TFS_REG, '<', $tfs_reg or die "Can't open TF motif region hit file, $tfs_reg: $!";
while (my $input = <$TFS_REG> ) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) [^\t]*? \t (\d+) /xms ) { 
        my $wbgene = $1;
        my $sites  = $2;
        # Because genes can have more than one regulatory sequence, it is crucial not to overwrite a high count with a lower count.
        if (     ( exists $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'}   )
             and ( $sites > $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'} ) ) {
            $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'} = $sites;
        }
        if (! exists $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'} ) {
            $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'} = $sites;
        }
    }
}
close $TFS_REG or die "Can't close filehandle to TF motif region hit file, $tfs_reg: $!";

open my $TFS_FLK, '<', $tfs_flk or die "Can't open TF motif flank hit file, $tfs_flk: $!";
while (my $input = <$TFS_FLK> ) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) [^\t]*? \t (\d+) /xms ) {
        my $wbgene = $1;
        my $sites  = $2;
        # Because genes can have more than one regulatory sequence, it is crucial not to overwrite a high count with a lower count.
        if (     ( exists $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'}   )
             and ( $sites > $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'} ) ) {
            $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'} = $sites;
        }
        if (! exists $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'} ) {
            $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'} = $sites;   
        }
    }
}
close $TFS_FLK or die "Can't close filehandle to TF motif flank hit file, $tfs_flk: $!";

open my $TFS_CONS_FLK, '<', $tfs_cons_flk or die "Can't open TF motif conserved flank hit file, $tfs_cons_flk: $!";
while (my $input = <$TFS_CONS_FLK> ) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) [^\t]*? \t (\d+) /xms ) {
        my $wbgene = $1;
        my $sites  = $2;
        # Because genes can have more than one regulatory sequence, it is crucial not to overwrite a high count with a lower count.
        if (     ( exists $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'}   ) 
             and ( $sites > $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'} ) ) {
            $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'} = $sites;
        }
        if (! exists $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'} ) { 
            $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'} = $sites;
        }
    }
}
close $TFS_CONS_FLK or die "Can't close filehandle to TF motif conserved flank hit file, $tfs_cons_flk: $!";

my $header = 'Gene' 
             . "\t" . 'Region_sites'
             . "\t" . 'Flank_sites'
             . "\t" . 'Cons_flank_sites'
             . "\t" . 'TF_type'
             . "\t" . 'Known_TreeFam'
             . "\t" . 'Possible_TreeFam'
             ;

my @wbgenes = sort keys %{ $data_ref->{'key_wbgene'} };
foreach my $wbgene (@wbgenes) { 
    my $region_sites = q{};
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'} ) {
        $region_sites = $data_ref->{'wbgene'}->{$wbgene}->{'region_sites'};
    }
    my $flank_sites = q{};
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'} ) {
        $flank_sites = $data_ref->{'wbgene'}->{$wbgene}->{'flank_sites'};
    }
    my $cons_flank_sites = q{};
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'} ) {
        $cons_flank_sites = $data_ref->{'wbgene'}->{$wbgene}->{'cons_flank_sites'};
    }
    my $tf_type = q{};
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'tf_type'} ) {
        $tf_type = $data_ref->{'wbgene'}->{$wbgene}->{'tf_type'};
    }
    my $known_treefam = q{};
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'known_treefam'} ) { 
        $known_treefam = $data_ref->{'wbgene'}->{$wbgene}->{'known_treefam'};
    }
    my $possible_treefam = q{};
    if ( exists $data_ref->{'wbgene'}->{$wbgene}->{'possible_treefam'} ) { 
        $possible_treefam = $data_ref->{'wbgene'}->{$wbgene}->{'possible_treefam'};
    }
    if ($header) { 
        print "$header\n";
        $header = q{};
    }
    print "$wbgene\t$region_sites\t$flank_sites\t$cons_flank_sites\t$tf_type\t$known_treefam\t$possible_treefam\n";
}

