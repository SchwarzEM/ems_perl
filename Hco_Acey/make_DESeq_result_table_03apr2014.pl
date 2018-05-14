#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $type = q{};

my $data_ref;

my @input_files = @ARGV;

foreach my $full_infile (@input_files) { 
    my $infile = basename $full_infile;
    # Sample: deseq_analyses/DESeq.comp.ACEY.L3i.vs.ACEY.24.PI.30mar2013.filt.txt or ...30mar2013.post_12.D_trip.filt.txt
    if ( $infile =~ /\A DESeq\.comp\.(ACEY\S+)\.vs\.(ACEY\S+)\.30mar2013(\S*)\.filt\.txt \z/xms ) { 
        my $cond_a = $1;
        my $cond_b = $2;
        $type      = $3;
        if ( ( $type ne q{} ) and ( $type ne '.post_12.D_trip' ) ) { 
             die "Can't parse type of file name $full_infile\n";
        }
        if ( $type eq q{} ) {
            $type = '(replicates: 17.D, 19.D)';
        }
        if ( $type eq '.post_12.D_trip' ) {
            $type = '(replicates: 12.D, 17.D, 19.D)';
        }
        my $transition = "$cond_a to $cond_b";
        open my $IN, '<', $full_infile;
        while (my $input = <$IN>) {
            chomp $input;
            # Sample inputs:
            # Acey_s0001.g113 7.98453841799826        3.23428293305673e-10
            # Acey_s0001.g115 -7.62781817670752       1.19779761223204e-06
            if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) { 
                my $gene  = $1;
                my $fold  = $2;
                my $q_val = $3;
                my $sign  = q{};

                if ($fold =~ /\A [-]/xms ) { 
                    $sign = 'down';
                }
                else {
                    $sign = 'up';
                }
                my $transition = "$cond_a to $cond_b [$sign]";
                if ( ( looks_like_number($q_val) ) and ( $q_val <= 0.01 ) ) { 
                    $data_ref->{'transition'}->{$transition}->{'gene'}->{$gene} = 1;
                }
            }
        }
    }
    else {
        die "Can't parse input file: $full_infile\n";
    }
}


my $header      = "Transition\tDESeq gene count $type";
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

