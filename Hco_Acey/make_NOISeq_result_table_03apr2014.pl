#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

my $data_ref;

my @input_files = @ARGV;

foreach my $full_infile (@input_files) { 
    my $infile = basename $full_infile;
    if ( $infile =~ /\A (AceyNoiseq\.(ACEY\S+)\.to\.(ACEY\S+)\.deg)\.(?:up|down)\.raw\.txt \z/xms ) { 
        my $stem   = $1;
        my $cond_a = $2;
        my $cond_b = $3;
        my $transition = "$cond_a to $cond_b";
        open my $IN, '<', $full_infile;
        while (my $input = <$IN>) {
            chomp $input;
            # Sample inputs:
            #   "ACEY.24.PI_mean" "ACEY.L3i_mean" "M" "D" "prob" "ranking"
            # "Acey_s0001.g130" 6252136.12558703 666.097749490496 13.1963276619642 6251470.02783754 1 6251470.02785147
            if ( ( $input !~ /ranking/ ) and ( $input =~ /\A \" (\S+) \" (?: \s+ \S+){2} \s+ (\S+) (?: \s+ \S+) \s+ \S+ \s+ \S+ \s* \z/xms ) ) { 
                my $gene  = $1;
                my $fold  = $2;
                my $sign  = q{};
                if ($fold > 0) {
                    $sign = 'up';
                }
                if ($fold < 0) {
                    $sign = 'down';
                }
                if ($fold == 0) {
                    die "Can't parse zero-fold log change of activity in input file $full_infile: $input\n";
                }
                my $transition = "$cond_a to $cond_b [$sign]";
                $data_ref->{'transition'}->{$transition}->{'gene'}->{$gene} = 1;
            }
        }
    }
    else {
        die "Can't parse input file: $full_infile\n";
    }
}

my $header      = "Transition\tNOISeq gene count";
my @transitions = sort keys %{ $data_ref->{'transition'} };

foreach my $transition (@transitions) { 
    my @genes      = sort keys %{ $data_ref->{'transition'}->{$transition}->{'gene'} };
    my $gene_count = @genes;
    $gene_count    = commify($gene_count);
    print "$header\n" if $header;
    $header = q{};
    print "$transition\t$gene_count\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

