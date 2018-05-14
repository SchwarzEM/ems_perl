#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $fasta = q{};
my $gene  = q{};
my $tx    = q{};

my $data_ref;

my $help;

GetOptions ( 'fasta=s' => \$fasta, 
             'help'    => \$help, );

if ($help or (! $fasta) ) { 
    die "Format: tabulate_ncoils_x.fa.pl --fasta|-f [ncoils-xed FASTA] --help|-h [print this message\n";
}

open my $FASTA, '<', $fasta or die "Can't open proteome: $fasta\n";
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ /\A > ((\S+\.g\d+)\.t\d+) /xms ) {
        $tx   = $1;
        $gene = $2;
    }
    elsif ( $tx and $gene and ( $input =~ /\S/xms ) ) { 
        $input =~ s/\s//g;
        $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx} .= $input;
    }
}
close $FASTA or die "Can't close filehandle to  proteome: $fasta\n";

my $header = "Gene\tNCoils\n";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @txs          = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'tx'} };
    my @coils_annots = ();
    foreach my $tx1 (@txs) { 
        my $sequence  = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1};
        my $len_seq   = length($sequence);
        my $len_coils = 0;
        if ( $sequence =~ /[x]+/xms ) { 
            $len_coils = ( $sequence =~ tr/x/x/ );
        }
        my $frac_coils = ($len_coils/$len_seq);
        $frac_coils    = sprintf "%.2f", $frac_coils;
        my $coils_text = "$frac_coils ($len_coils/$len_seq)";
        if ( $frac_coils > 0 ) { 
            push @coils_annots, $coils_text;
        }
    }
    my $coils_annot = join '; ', @coils_annots;
    print $header if $header;
    $header = q{};
    if ( $coils_annot =~ /\S/xms ) { 
        $coils_annot = 'NCoils: ' . $coils_annot;
    }
    print "$gene1\t$coils_annot\n";
}
