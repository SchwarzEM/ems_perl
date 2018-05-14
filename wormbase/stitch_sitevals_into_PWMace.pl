#!/usr/bin/env perl

# stitch_sitevals_into_PWMace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/25/2008.
# Purpose: move site values from raw PWM.ace (with good ones) to mature PWM.ace (with bad ones).

use strict;
use warnings;

my $old_obj2_textA_ref;
my $old_obj2_textB_ref;

my %old_obj2new_obj;
my $new_obj2_sites_ref;

my $old_object = q{};
my $new_object = q{};
my $read_textB = 0;

while (my $input = <>) { 
    if ( $input =~ / \A 
                     Position_Matrix 
                     \s+ : \s+ 
                     \"
                     (WBPmat\d+)\.pwm
                     \"
                   /xms ) { 
        $new_object = $1;
        $old_object = q{};
        $read_textB = 0;
    }
    elsif ( $input =~ / \A
                        Position_Matrix
                        \s+ : \s+
                        \"
                        (WBPmat\d+)
                        \"
                     /xms ) { 
        $old_object = $1;
        $new_object = q{};
        $read_textB = 0;
        push @{ $old_obj2_textA_ref->{$old_object} }, $input;
    }
    elsif (     ( $old_object       ) 
            and ( $read_textB == 0  ) 
            and ( $input =~ /\S/xms ) 
            and ( $input !~ /\A \s* Site_values /xms ) ) { 
        push @{ $old_obj2_textA_ref->{$old_object} }, $input;
    }
    elsif (     ( $old_object       )
            and (! $read_textB      )   
            and ( $input =~ /\A \s* Site_values /xms ) ) {
        $read_textB = 1;
    }
    elsif (     ( $old_object       )
            and ( $read_textB       )
            and ( $input =~ /\S/xms )
            and ( $input !~ /\A \s* Site_values /xms ) ) {
        push @{ $old_obj2_textB_ref->{$old_object} }, $input;
        if ( $input =~ / \A
                         Remark .+
                         from \s+ frequency \s+ matrix \s+
                         (WBPmat\d+)
                       /xms ) { 
            my $new = $1;
            $old_obj2new_obj{$old_object} = $new;
        }
    }
    elsif ( ( $new_object ) and ( $input =~ /\A \s* Site_values /xms ) ) {
        push @{ $new_obj2_sites_ref->{$new_object} }, $input;
    }
}

print "\n";
foreach my $obj_w_most_stuff (sort keys %{ $old_obj2_textA_ref } ) { 
    if ( exists $old_obj2new_obj{$obj_w_most_stuff} ) { 
        my $obj_w_site_values = $old_obj2new_obj{$obj_w_most_stuff};
        foreach my $output1 ( @{ $old_obj2_textA_ref->{$obj_w_most_stuff} } ) { 
            print $output1;
        }
        foreach my $output2 ( @{ $new_obj2_sites_ref->{$obj_w_site_values} } ) { 
            print $output2;
        }
        if ( exists $old_obj2_textB_ref->{$obj_w_most_stuff} ) { 
            foreach my $output3 ( @{ $old_obj2_textB_ref->{$obj_w_most_stuff} } ) { 
                print $output3;
            }
        }
        print "\n";
    }
}


