#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $orig_table = q{};
my $sec_list   = q{};

$orig_table = $ARGV[0] if $ARGV[0];
$sec_list   = $ARGV[1] if $ARGV[1];

my %sec_genes = ();

if ( (! $orig_table) or (! $sec_list) ) {
    die "Format: add_secretome_annots_15oct2018.pl [orig annot table] [secretome genelist] > [edited annot table]\n";
}

open my	$SEC, '<', $sec_list;
while (my $input = <$SEC>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        $sec_genes{$input} = 1;
    }
    else {
        die "From secretome gene list $sec_list, cannot parse: $input\n";
    }
}
close $SEC;

open my $ORIG, '<', $orig_table;
while (my $input = <$ORIG>) {
    chomp $input;
    if ( $input =~ /\A ( (\S+) (?: \t [^\t]*){4} \t ) ([^\t]*) (\t .*) \z/xms ) {
        my $text1   = $1;
        my $gene    = $2;
        my $phobius = $3;
        my $text2   = $4;
        if ( exists $sec_genes{$gene} ) {
            if ( $phobius =~ /\S/xms ) { 
                $phobius = "$phobius [+SecP]";
                $phobius =~ s/[ ]+/ /g;
            }
            else {
                $phobius =  'SecretomeP';
            }
        }
        print "$text1$phobius$text2\n";
    }
    else {
        die "From original annot. table $orig_table, cannot parse: $input\n";
    }
}
close $ORIG;

