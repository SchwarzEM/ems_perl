#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $seg      = q{};
my $proteome = q{};

my $gene = q{};
my $tx   = q{};

my $data_ref;
my $help;

GetOptions ( 'seg=s'      => \$seg,
             'proteome=s' => \$proteome,
             'help'       => \$help, );

if ( $help or (! $seg) or (! $proteome) ) { 
    die "Format: summarize_psegs_08nov2012.pl --seg|-s [FASTA of low-complexity pseg domains] --proteome|-p [proteome; both -s and -p assume AUGUSTUS gene/trx names]\n";
}

open my $PROTEOME, '<', $proteome or die "Can't open proteome: $proteome\n";
while (my $input = <$PROTEOME>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ((\S+\.g\d+)\.t\d+) /xms ) { 
            $tx   = $1;
            $gene = $2;
        }
        else { 
            die "Can't parse: $input\n";
        }
    }
    else {
        $input =~ s/\s//g;
        my $len = length($input);
        $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'len'} += $len;
    }
}
close $PROTEOME or die "Can't close filehandle to proteome: $proteome\n";

open my $PSEG, '<', $seg or die "Can't open filehandle to PSEG FASTA of low-complexity domains: $seg\n";
while (my $input = <$PSEG>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ((\S+\.g\d+)\.t\d+) \( (\d+) \- (\d+) \) /xms ) { 
            $tx       = $1;
            $gene     = $2;
            my $start = $3;
            my $stop  = $4;
            my $seg_span = $stop - $start + 1;
            $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'seg_span'} += $seg_span;
        }
        else { 
            die "Can't parse: $input\n";
        }
    }
}
close $PSEG or die "Can't close filehandle to PSEG FASTA of low-complexity domains: $seg\n";

my @genes = grep { /\S/ } sort keys %{ $data_ref->{'gene'} };
my $header = "Gene\tPsegs\n";
foreach my $gene1 (@genes) { 
    my @annots = q();
    my @txs    = grep { /\S/ } sort keys %{ $data_ref->{'gene'}->{$gene1}->{'tx'} };
    foreach my $tx1 (@txs) { 
        if ( exists $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'len'}) { 
            my $length   = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'len'};
            my $seg_span = 0;
            if ( exists $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'seg_span'}) { 
                $seg_span = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'seg_span'};
            }
            if ( $seg_span > 0 ) { 
                my $fraction = ($seg_span/$length);
                $fraction = sprintf "%.2f", $fraction;
                my $annot = "$fraction ($seg_span/$length)";
                push @annots, $annot;
            }
        }
        else { 
            die "Failed to get length for protein product of $tx1\n";
        }
    }
    my $annot_text = q{};
    if (@annots) { 
        @annots = grep { /\S/ } @annots;
        $annot_text = join '; ', @annots;
    }
    print $header if $header;
    $header = q{};
    if ( $annot_text =~ /\S/xms ) { 
        $annot_text = 'Seg: ' . $annot_text;
    }
    print "$gene1\t$annot_text\n";
}

