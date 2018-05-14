#!/usr/bin/env perl

# tabulate_phobius_hits.pl -- Erich Schwarz <ems394@cornell.edu>, 9/28/2014.
# Given a gene-to-CDS table and an CDS-based Phobius short output, produce a *.tsv.txt table linking genes to summarized Phobius annotations.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw{uniq};
use autodie;

my $gene2cds = q{};
my $phobius  = q{};

my $data_ref;

my $help;

GetOptions ( 'gene2cds=s' => \$gene2cds,
             'phobius=s'  => \$phobius, 
             'help'       => \$help, );

if ($help or (! $gene2cds ) or (! $phobius) ) { 
    die "Format: tabulate_phobius_hits.pl\n",
        "    --gene2cds|-g  [gene to CDS name table, tab-delimited]\n",
        "    --phobius|-p   [phobius short-output file; has CDS names and SP/TM results]\n",
        "    --help|-h      [prints this message]\n",
        ;
}

open my $GENE, '<', $gene2cds;
while (my $input = <$GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gene = $1;
        my $cds  = $2;
        if ( ( exists $data_ref->{'cds'}->{$cds}->{'gene'} ) and ( $data_ref->{'cds'}->{$cds}->{'gene'} ne $gene ) ) {
            die "CDS $cds is inconsistently mapped to both $data_ref->{'cds'}->{$cds}->{'gene'} and $gene\n";
        }
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else {
        die "From gene-to-CDS name table $gene2cds, can't parse input line: $input\n";
    }
}
close $GENE;

open my $PHOB, '<', $phobius;
while (my $input = <$PHOB>) { 
    chomp $input;
    if ( $input !~ /\A SEQENCE \s+ ID \s+ TM \s+ SP \s + PREDICTION/xms ) { 

        # Typical input:
        # 
        # Acey_2012.08.05_0001.g10.t2     1  Y n2-13c18/19o136-153i
        # Acey_2012.08.05_0001.g10.t3     0  Y n2-13c18/19o
        # Acey_2012.08.05_0001.g100.t1    0  0 o

        if ( $input =~ /\A (\S+) \s+ (\d+) \s+ (0|Y) \s+ \S+ \z /xms ) { 
            my $cds       = $1;
            my $tmem_no = $2;
            my $sig_seq = $3;
            my @phob_annots = ();
            my $phob_annot_text = q{};
            if ( $sig_seq eq 'Y' ) { 
                push @phob_annots, 'SigP';
            }
            if ( $tmem_no >= 1 ) { 
                my $tm_text = 'TM(' . $tmem_no . 'x)';
                push @phob_annots, $tm_text;
            }
            if (@phob_annots) {
                $phob_annot_text = join '+', @phob_annots;
            }
            if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
                die "Can't link CDS $cds to gene\n";
            }
            my $gene = $data_ref->{'cds'}->{$cds}->{'gene'};
            $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds}->{'phob_annot'} = $phob_annot_text;
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
}
close $PHOB;

my $header = "Gene\tPhobius\n";
my @genes  = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @phob_annots = ();
    my @cdses = sort keys %{ $data_ref->{'gene'}->{$gene}->{'cds'} };
    foreach my $cds (@cdses) {
        my $init_phob_annot_text = q{};
        if ( $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds}->{'phob_annot'} !~ /\A\s*\z/xms ) {
            $init_phob_annot_text = $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds}->{'phob_annot'};
            push @phob_annots, $init_phob_annot_text;
        }
    }
    if (@phob_annots) { 
        @phob_annots = sort @phob_annots;
        @phob_annots = uniq @phob_annots;
        my $phob_annot_text = join '; ', @phob_annots;
        print $header if $header;
        $header = q{};
        print "$gene\t$phob_annot_text\n";
    }
}

