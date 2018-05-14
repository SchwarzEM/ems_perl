#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $acey_data = $ARGV[0];
my $nam_data  = $ARGV[1];
my $cel_data  = $ARGV[2];

my $data_ref;

my %ok_repeat_types = ( 'DF0000212.2' => 1,  # HSMAR1
                        'DF0000213.2' => 1,  # HSMAR2
                        'DF0000362.2' => 1,  # L3
                        'DF0001135.2' => 1,  # Plat_L3
);

my $header = "Element_type\tElement_accession\tE-value\tSpecies\tRepeat_name\tElement_description";

open my $ACEY, '<', $acey_data;
while (my $input = <$ACEY>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ (?: \S+ \s+){9} (\S+) \s+ (?: \S+ \s+){2} (.+) \z/xms ) { 
        my $target_name = $1;
        my $accession   = $2;
        my $query_name  = $3;
        my $e_value     = $4;
        my $description = $5;
        if ( exists $ok_repeat_types{$accession} ) { 
            my $out_text = "$target_name\t$accession\t$e_value\tA. ceylanicum\t$query_name\t$description";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "A. ceylanicum\t$query_name";
        }
    }
}
close $ACEY;

open my $NAM, '<', $nam_data;
while (my $input = <$NAM>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ (?: \S+ \s+){9} (\S+) \s+ (?: \S+ \s+){2} (.+) \z/xms ) {
        my $target_name = $1;
        my $accession   = $2;
        my $query_name  = $3;
        my $e_value     = $4;
        my $description = $5;
        if ( exists $ok_repeat_types{$accession} ) {
            my $out_text = "$target_name\t$accession\t$e_value\tN. americanus\t$query_name\t$description";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "N. americanus\t$query_name";
        }
    }
}
close $NAM;

open my $CEL, '<', $cel_data;
while (my $input = <$CEL>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ (?: \S+ \s+){9} (\S+) \s+ (?: \S+ \s+){2} (.+) \z/xms ) {
        my $target_name = $1;
        my $accession   = $2;
        my $query_name  = $3;
        my $e_value     = $4;
        my $description = $5;    
        if ( exists $ok_repeat_types{$accession} ) {
            my $out_text = "$target_name\t$accession\t$e_value\tC. elegans\t$query_name\t$description";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "C. elegans\t$query_name";
        }
    }
}
close $CEL;

my @out_lines = sort { $data_ref->{'out_text'}->{$a}->{'e_value'} <=> $data_ref->{'out_text'}->{$b}->{'e_value'} } keys %{ $data_ref->{'out_text'} };

print "$header\n";

foreach my $out_line (@out_lines) {
    my $element = $data_ref->{'out_text'}->{$out_line}->{'element'};
    if (! exists $data_ref->{'printed'}->{$element} ) { 
        print "$out_line\n";
    }
    $data_ref->{'printed'}->{$element} = 1;
}


