#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $phobius = q{};
my $gene    = q{};
my $tx      = q{};

my $data_ref;

my $help;

GetOptions ( 'phobius=s' => \$phobius, 
             'help'      => \$help, );

if ($help or (! $phobius) ) { 
    die "Format: tabulate_phobius_hits_14nov2012.pl\n",
        "    --phobius|-p  [phobius short-output file; has gene/trx names and SP/TM results]\n",
        "    --help|-h     [prints this message]\n",
        ;
}

open my $PHOB, '<', $phobius or die "Can't open proteome: $phobius\n";
while (my $input = <$PHOB>) { 
    chomp $input;
    if ( $input !~ /\A SEQENCE \s+ ID \s+ TM \s+ SP \s + PREDICTION/xms ) { 

        # Typical input:
        # 
        # Acey_2012.08.05_0001.g10.t2     1  Y n2-13c18/19o136-153i
        # Acey_2012.08.05_0001.g10.t3     0  Y n2-13c18/19o
        # Acey_2012.08.05_0001.g100.t1    0  0 o

        if ( $input =~ /\A ((\S+\.g\d+)\.t\d+) \s+ (\d+) \s+ (0|Y) \s+ \S+ \z /xms ) { 
            $tx         = $1;
            $gene       = $2;
            my $tmem_no = $3;
            my $sig_seq = $4;
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
            $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'phob_annot'} = $phob_annot_text;
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
}
close $PHOB or die "Can't close filehandle to proteome: $phobius\n";

my $header = "Gene\tPhobius\n";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @txs         = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'tx'} };
    my %phob_annots = ();
    my $phob_annot_text = q{};
    foreach my $tx1 (@txs) {
        if ( $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'phob_annot'} !~ /\A\s*\z/xms ) { 
            my $init_phob_annot = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'phob_annot'};
            $phob_annots{$init_phob_annot} = 1;
        }
    }
    if (%phob_annots) { 
        my @phob_annot_list = sort keys %phob_annots;
        $phob_annot_text = join '; ', @phob_annot_list;
    }
    print $header if $header;
    $header = q{};
    print "$gene1\t$phob_annot_text\n";
}

