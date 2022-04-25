#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $data_ref;

my $infile = $ARGV[0];
my $thresh = $ARGV[1];

my $header = "Type\tProt.-cod.\tncRNA";

if ( (! $infile) or (! $thresh ) ) {
    die "Format: gene_expr_thresh_2022.04.25.01.pl [infile] [expr. threshold] > [genes per condition]\n";
}

if (! looks_like_number($thresh) ) {
    warn "Format: gene_expr_thresh_2022.04.25.01.pl [infile] [expr. threshold] > [genes per condition]\n";
    die "[expr. threshold] must be a number!\n";
}

my $type_count = 0;

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;

    if ( $input =~ /\A [^\t\s]+ \t ([^\t]+) \t (.+) \z/xms ) {
        my $gene_annot = $1;
        my $type_text  = $2;

        if ( $gene_annot =~ /protein/xms ) {
            $gene_annot = 'Protein-coding';
        }

        my @types = split /\t/, $type_text;
        $type_count = @types;
        my $array_index = $type_count;
        $array_index--;
        foreach my $i (0..$array_index) {
            my $j = ($i + 1);
            push @{ $data_ref->{$j}->{$gene_annot} }, $types[$i];
        }
    }
}

foreach my $i (1..$type_count) {

    my @protein_data = ();
    my @ncRNA_data = ();
    my $type = q{};

    if ( exists $data_ref->{$i}->{'Protein-coding'} ) {
        @protein_data = @{ $data_ref->{$i}->{'Protein-coding'} };
    }

    if ( exists $data_ref->{$i}->{'ncRNA'} ) {
        @ncRNA_data = @{ $data_ref->{$i}->{'ncRNA'} };
    }

    if ( exists $data_ref->{$i}->{'Coding'} ) {
        $type = shift @{ $data_ref->{$i}->{'Coding'} };
    }

    my @protein_expr = grep { $_ >= $thresh } @protein_data;
    my $protein_expr_count = @protein_expr;
    $protein_expr_count = commify($protein_expr_count);

    my @ncRNA_expr = grep { $_ >= $thresh } @ncRNA_data;
    my $ncRNA_expr_count = @ncRNA_expr;
    $ncRNA_expr_count = commify($ncRNA_expr_count);

    print "$header\n" if $header;
    $header = q{};

    print "$type\t$protein_expr_count\t$ncRNA_expr_count\n";
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
