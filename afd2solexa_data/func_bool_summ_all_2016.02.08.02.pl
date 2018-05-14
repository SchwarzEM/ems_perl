#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my @input_files = @ARGV;

my $data_ref;

my $header = "GO_term\tp-value[condition]";

foreach my $infile (@input_files) {
    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) {
        chomp $input;
        # sample input line:
        # defense response to bacterium [GO:0042742]	9.05405e-12
        if ( $input =~ /\A (\S.*\S) \s+ \[ (GO:\d+) \] \t (\S+) \z/xms ) {
            my $go_desc = $1;
            my $go_term = $2;
            my $p_value = $3;

            if (! looks_like_number($p_value) ) {
                die "p-value appears non-numerical in: $input\n";
            }
            if (     ( exists $data_ref->{'go_term'}->{$go_term}->{'go_desc'}      ) 
                 and ( $data_ref->{'go_term'}->{$go_term}->{'go_desc'} ne $go_desc ) ) {
                die "GO term $go_term has two different descriptions: $data_ref->{'go_term'}->{$go_term}->{'go_desc'} vs. $go_desc\n";
            }

            $data_ref->{'go_term'}->{$go_term}->{'go_desc'} = $go_desc;
            $data_ref->{'go_term'}->{$go_term}->{'p_value'}->{$p_value}->{'condition'}->{$infile} = 1;
            $data_ref->{'p_value'}->{$p_value}->{'go_term'}->{$go_term} = 1;
        }
        elsif ( $input !~ /\A GO_term \t p-value \z/xms ) { 
            die "In input file $infile, cannot parse: $input\n";
        }
    }
    close $INFILE;
}

my @p_values = sort { $a <=> $b } keys %{ $data_ref->{'p_value'} };
foreach my $p_value_best (@p_values) {
     my @go_terms = sort keys %{ $data_ref->{'p_value'}->{$p_value_best}->{'go_term'} };
     foreach my $go_term (@go_terms) {
         my $go_desc      = $data_ref->{'go_term'}->{$go_term}->{'go_desc'};
         my $go_text      = "$go_desc [$go_term]";

         my @p_values_all = sort { $a <=> $b } keys %{ $data_ref->{'go_term'}->{$go_term}->{'p_value'} };

         my @p_conds_all  = ();
         foreach my $p_value_all (@p_values_all) {
             my @conditions_all = sort keys %{ $data_ref->{'go_term'}->{$go_term}->{'p_value'}->{$p_value_all}->{'condition'} };         

             my $condition_text = join ', ', @conditions_all;
             my $p_cond_all     = "$p_value_all [$condition_text]";
             push @p_conds_all, $p_cond_all;

             # Once I have used these data, delete them, so that that I will only access them once:
             delete $data_ref->{'go_term'}->{$go_term}->{'p_value'}->{$p_value_all}->{'condition'};
         }
         my $p_cond_text = join '; ', @p_conds_all;

         # Prevent empty reprints of GO terms by requiring there to be some remaining data.
         if ( $p_cond_text =~ /\S/xms ) { 
             # Print header exactly once.
             print "$header\n" if $header;
             $header = q{};

             print "$go_text\t$p_cond_text\n";
         }

         # Delete the data once I have used them:
         delete $data_ref->{'go_term'}->{$go_term}->{'p_value'};
     }

     # Delete the data once I have used them:
     delete $data_ref->{'p_value'}->{$p_value_best}->{'go_term'};
}

