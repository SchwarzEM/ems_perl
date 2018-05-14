#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $acey_data = $ARGV[0];
my $nam_data  = $ARGV[1];
my $cel_data  = $ARGV[2];
my $cbri_data = $ARGV[3];

my $data_ref;

my %ok_repeat_types = ( 'DF0000110.2' => 1,  # CR1_Mam
                        'DF0000212.2' => 1,  # HSMAR1
                        'DF0000213.2' => 1,  # HSMAR2
                        'DF0000362.2' => 1,  # L3
                        'DF0001135.2' => 1,  # Plat_L3
);

my $header = "Element_type\tElement_accession\tE-value\tSpecies\tRepeat_name";

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
            my $out_text = "$target_name\t$accession\t$e_value\tA. ceylanicum\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "A. ceylanicum\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'accession'} = $accession;
            # Get the lowest E-value observed for each accession:
            if (    (! exists $data_ref->{'accession'}->{$accession}->{'e_value'} ) 
                 or ( $data_ref->{'accession'}->{$accession}->{'e_value'} > $e_value ) ) {
                $data_ref->{'accession'}->{$accession}->{'e_value'} = $e_value;
            }
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
            my $out_text = "$target_name\t$accession\t$e_value\tN. americanus\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "N. americanus\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'accession'} = $accession;
            # Get the lowest E-value observed for each accession:
            if (    (! exists $data_ref->{'accession'}->{$accession}->{'e_value'} )
                 or ( $data_ref->{'accession'}->{$accession}->{'e_value'} > $e_value ) ) {
                $data_ref->{'accession'}->{$accession}->{'e_value'} = $e_value;
            }
        }
    }
}
close $NAM or die;

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
            my $out_text = "$target_name\t$accession\t$e_value\tC. elegans\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "C. elegans\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'accession'} = $accession;
            # Get the lowest E-value observed for each accession:
            if (    (! exists $data_ref->{'accession'}->{$accession}->{'e_value'} )
                 or ( $data_ref->{'accession'}->{$accession}->{'e_value'} > $e_value ) ) {
                $data_ref->{'accession'}->{$accession}->{'e_value'} = $e_value;
            }
        }
    }
}
close $CEL or die;

open my $CBRI, '<', $cbri_data;
while (my $input = <$CBRI>) {
    chomp $input; 
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ (?: \S+ \s+){9} (\S+) \s+ (?: \S+ \s+){2} (.+) \z/xms ) {
        my $target_name = $1;
        my $accession   = $2;
        my $query_name  = $3;
        my $e_value     = $4;
        my $description = $5;
        if ( exists $ok_repeat_types{$accession} ) {  
            my $out_text = "$target_name\t$accession\t$e_value\tC. briggsae\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'e_value'} = $e_value;
            $data_ref->{'out_text'}->{$out_text}->{'element'} = "C. briggsae\t$query_name";
            $data_ref->{'out_text'}->{$out_text}->{'accession'} = $accession;
            # Get the lowest E-value observed for each accession:
            if (    (! exists $data_ref->{'accession'}->{$accession}->{'e_value'} )
                 or ( $data_ref->{'accession'}->{$accession}->{'e_value'} > $e_value ) ) {
                $data_ref->{'accession'}->{$accession}->{'e_value'} = $e_value;
            }
        }
    }
}
close $CBRI or die;

my @out_lines_initial = sort { $data_ref->{'out_text'}->{$a}->{'e_value'} <=> $data_ref->{'out_text'}->{$b}->{'e_value'} } keys %{ $data_ref->{'out_text'} };

# Censor all out lines that have suboptimal E-values:
foreach my $out_line_initial (@out_lines_initial) {
    my $element = $data_ref->{'out_text'}->{$out_line_initial}->{'element'};
    if ( exists $data_ref->{'printed'}->{$element} ) {
        delete $data_ref->{'out_text'}->{$out_line_initial};
    }
    else { 
        $data_ref->{'printed'}->{$element} = 1;
    }
}

# Having done such censorship across the entire list as a whole, ordered by ascending E-values, it is now safe to split this into lists that 
#    are more intelligible because they are presented by *element* type.

my @out_lines = sort { $data_ref->{'out_text'}->{$a}->{'e_value'} 
                       <=> $data_ref->{'out_text'}->{$b}->{'e_value'} } 
                keys %{ $data_ref->{'out_text'} };

my @accessions_to_list = sort { $data_ref->{'accession'}->{$a}->{'e_value'} 
                                <=> $data_ref->{'accession'}->{$b}->{'e_value'} } 
                         keys %{ $data_ref->{'accession'} };

foreach my $acc_to_list (@accessions_to_list) { 
    print "$header\n";
    foreach my $out_line (@out_lines) {
        my $element = $data_ref->{'out_text'}->{$out_line}->{'element'};
        my $accession = $data_ref->{'out_text'}->{$out_line}->{'accession'};
        if ( $accession eq $acc_to_list ) { 
            print "$out_line\n";
        }
    }
    print "\n";
}


