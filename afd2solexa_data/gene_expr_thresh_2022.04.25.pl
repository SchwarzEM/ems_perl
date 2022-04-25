#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $data_ref;

my $infile = $ARGV[0];
my $thresh = $ARGV[1];

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
    if ( $input =~ /\A [^\t\s]+ \t (?: Coding |[^\t]* protein [^\t]*) \t (.+) \z/xms ) {
        my $type_text = $1;
        my @types = split /\t/, $type_text;
        $type_count = @types;
        my $array_index = $type_count;
        $array_index--;
        foreach my $i (0..$array_index) {
            my $j = ($i + 1);
            push @{ $data_ref->{$j} }, $types[$i];
        }
    }
}

foreach my $i (1..$type_count) {
    my @data = @{ $data_ref->{$i} };
    my $type = shift @data;
    my @expr = grep { $_ >= $thresh } @data;
    my $expr_count = @expr;
    $expr_count = commify($expr_count);
    print "$type\t$expr_count\n";
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
