#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my $data_ref;

my @transitions = (
    'L3i_to_24.PI',
    'L3i_to_24.HCM',
    '24.PI_to_5.D',
    '5.D_to_12.D',
    '12.D_to_17.D',
    '17.D_to_19.D',
);

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A ([^\t]+ \t) [^\t]* \t ([^\t]+ \t [^\t]+) \t ([^\t]+) \z/xms ) {   
        my $go_term1 = $1;
        my $go_term2 = $2;
        my $go_stats = $3;
        my $go_term  = $go_term1 . $go_term2;

        $data_ref->{'GO_term'}->{$go_term}->{'L3i_to_24.PI'}  = 1;
        $data_ref->{'GO_term'}->{$go_term}->{'L3i_to_24.HCM'} = 1;
        $data_ref->{'GO_term'}->{$go_term}->{'24.PI_to_5.D'}  = 1;
        $data_ref->{'GO_term'}->{$go_term}->{'5.D_to_12.D'}   = 1;
        $data_ref->{'GO_term'}->{$go_term}->{'12.D_to_17.D'}  = 1;
        $data_ref->{'GO_term'}->{$go_term}->{'17.D_to_19.D'}  = 1;

        if ( $go_stats =~ / Acey_L3i [ ] to [ ] Acey_24.PI [ ] \[ [+] : [ ] (\S+) \] ;/xms ) { 
            my $p_value = $1;
            if (! looks_like_number($p_value) ) {
                die "Can't parse putative p-value: $p_value\n";
            }
            $data_ref->{'OK_GO_term'}->{$go_term} = 1;
            $data_ref->{'GO_term'}->{$go_term}->{'L3i_to_24.PI'} = $p_value;
        }
        # Cope with 'Acey_24HCM' typo, rather than the 'Acey_24.HCM' which would have been preferable and consistent nomenclature:
        if ( $go_stats =~ / Acey_L3i [ ] to [ ] Acey_24HCM [ ] \[ [+] : [ ] (\S+) \] ;/xms ) {
            my $p_value = $1;
            if (! looks_like_number($p_value) ) {
                die "Can't parse putative p-value: $p_value\n";
            }
            $data_ref->{'OK_GO_term'}->{$go_term} = 1;
            $data_ref->{'GO_term'}->{$go_term}->{'L3i_to_24.HCM'} = $p_value;
        }
        if ( $go_stats =~ / Acey_24.PI [ ] to [ ] Acey_5.D [ ] \[ [+] : [ ] (\S+) \] ;/xms ) { 
            my $p_value = $1;
            if (! looks_like_number($p_value) ) {
                die "Can't parse putative p-value: $p_value\n";
            }
            $data_ref->{'OK_GO_term'}->{$go_term} = 1;
            $data_ref->{'GO_term'}->{$go_term}->{'24.PI_to_5.D'} = $p_value;
        }
        if ( $go_stats =~ / Acey_5.D [ ] to [ ] Acey_12.D [ ] \[ [+] : [ ] (\S+) \] ;/xms ) {
            my $p_value = $1;
            if (! looks_like_number($p_value) ) {
                die "Can't parse putative p-value: $p_value\n";
            }
            $data_ref->{'OK_GO_term'}->{$go_term} = 1;
            $data_ref->{'GO_term'}->{$go_term}->{'5.D_to_12.D'} = $p_value;
        }
        if ( $go_stats =~ / Acey_12.D [ ] to [ ] Acey_17.D [ ] \[ [+] : [ ] (\S+) \] ;/xms ) {
            my $p_value = $1;
            if (! looks_like_number($p_value) ) {
                die "Can't parse putative p-value: $p_value\n";
            }
            $data_ref->{'OK_GO_term'}->{$go_term} = 1;
            $data_ref->{'GO_term'}->{$go_term}->{'12.D_to_17.D'} = $p_value;
        }
        if ( $go_stats =~ / Acey_17.D [ ] to [ ] Acey_19.D [ ] \[ [+] : [ ] (\S+) \] ;/xms ) {
            my $p_value = $1; 
            if (! looks_like_number($p_value) ) {
                die "Can't parse putative p-value: $p_value\n";
            }
            $data_ref->{'OK_GO_term'}->{$go_term} = 1;
            $data_ref->{'GO_term'}->{$go_term}->{'17.D_to_19.D'} = $p_value;
        }
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

# The only GO terms we care about are the ones that had some significance:
my @ok_go_terms = sort keys %{ $data_ref->{'OK_GO_term'}  };

foreach my $go_term (@ok_go_terms) { 
    foreach my $transition (@transitions) {
        push @{ $data_ref->{'GO_term'}->{$go_term}->{'p_values'} }, $data_ref->{'GO_term'}->{$go_term}->{$transition};
    }
    my @p_values = @{ $data_ref->{'GO_term'}->{$go_term}->{'p_values'} };

    # Resulting sort is like this: 7.84166e-10 1.39678e-06 1 1 1;
    @p_values = sort { $a <=> $b } @p_values;

    # ... then parse their meaning into p-value-sortable version.
    foreach my $transition (@transitions) {
        if ( $data_ref->{'GO_term'}->{$go_term}->{$transition} == $p_values[0] ) { 
            my $p_val_output = $data_ref->{'GO_term'}->{$go_term}->{$transition};
            my $output_text  = "$transition\t$go_term\t$p_val_output";
            $data_ref->{'p_val_output'}->{$p_val_output}->{'output_text'}->{$output_text} = 1;
        }
    }
}

my @p_val_outputs = sort { $a <=> $b } keys %{ $data_ref->{'p_val_output'} };

foreach my $p_val_output (@p_val_outputs) {
    my @output_texts = sort keys %{ $data_ref->{'p_val_output'}->{$p_val_output}->{'output_text'} };
    foreach my $output_text (@output_texts) {
        print "$output_text\n";
    }
}

