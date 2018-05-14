#!/usr/bin/env perl

use strict;
use warnings;

my $acey_data = $ARGV[0];
my $cel_data  = $ARGV[1];

my %ok_repeat_types = ( 'DF0000212.2' => 1,  # HSMAR1
                        'DF0000213.2' => 1,  # HSMAR2
                        'DF0001135.2' => 1,  # Plat_L3
);

my %output_lines = ();

my $header = "Element_type\tElement_accession\tE-value\tSpecies\tRepeat_name";


open my $ACEY, '<', $acey_data or die "Can't open A. ceylanicum DFAM data $acey_data: $!";
while (my $input = <$ACEY>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ (?: \S+ \s+){9} (\S+) \s+ (?: \S+ \s+){2} (.+) \z/xms ) { 
        my $target_name = $1;
        my $accession   = $2;
        my $query_name  = $3;
        my $e_value     = $4;
        my $description = $5;
        if ( exists $ok_repeat_types{$accession} ) { 
            my $out_text = "$target_name\t$accession\t$e_value\tA. ceylanicum\t$query_name";
            $output_lines{$out_text} = $e_value;
        }
    }
}
close $ACEY or die "Can't close filehandle to A. ceylanicum DFAM data $acey_data: $!";

open my $CEL, '<', $cel_data or die "Can't open C. elegans DFAM data $cel_data: $!";
while (my $input = <$CEL>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ (?: \S+ \s+){9} (\S+) \s+ (?: \S+ \s+){2} (.+) \z/xms ) {
        my $target_name = $1;
        my $accession   = $2;
        my $query_name  = $3;
        my $e_value     = $4;
        my $description = $5;    
        if ( exists $ok_repeat_types{$accession} ) {
            my $out_text = "$target_name\t$accession\t$e_value\tC. elegans\t$query_name";
            $output_lines{$out_text} = $e_value;
        }
    }
}
close $CEL or die "Can't close filehandle to C. elegans DFAM data $cel_data: $!";

my @out_lines = sort { $output_lines{$a} <=> $output_lines{$b} } keys %output_lines;

print "$header\n";

foreach my $out_line (@out_lines) {
    print "$out_line\n";
}

