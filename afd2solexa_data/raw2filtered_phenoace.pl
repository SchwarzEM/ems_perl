#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my $var_rnai      = q{};
my $curr_var_rnai = q{};
my $gene          = q{};
my $pheno         = q{};
my $conf          = q{};
my $evidence      = q{};
my $print_start   = 0;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A (?:RNAi|Variation) \s+ : \s+ \" ([^\"]+) \" /xms ) { 
        $var_rnai = $1;
        $data_ref->{$var_rnai}->{'header'} = $input;

        # Big, annoying analysis and printout subroutine:
        export_and_clear_stored_data();

        # Avoid inadvertant variable carry-over:
        $curr_var_rnai = q{};
        $gene =          q{};
    }
    elsif ( $input =~ / \A Gene \s+ \" (WBGene\d+) \" /xms ) { 
        $gene = $1;
        $curr_var_rnai = $var_rnai;
        $data_ref->{$curr_var_rnai}->{'gene'}->{$gene} = $input;
        $var_rnai = q{};
    }
    # Phenotype        "WBPhenotype:0000062" Curator_confirmed "WBPerson712"
    # Phenotype        "WBPhenotype:0000062" Not Curator_confirmed "WBPerson712"
    elsif ( $input =~ /\A Phenotype \s+ \" WBPhenotype: (\d+) \" \s+ Curator_confirmed \s+ \" (WBPerson\d+) \" /xms ) {
        $pheno = $1;
        $conf  = $2;
        $data_ref->{$curr_var_rnai}->{'pgene'}->{$gene}->{'pheno'}->{$pheno}->{'conf'}->{$conf}->{'yes'} = $input;
    }
    elsif ( $input =~ /\A Phenotype \s+ \" WBPhenotype: (\d+) \" \s+ Not \s+ Curator_confirmed \s+ \" (WBPerson\d+) \" /xms ) {
        $pheno = $1;
        $conf  = $2;
        $data_ref->{$curr_var_rnai}->{'pgene'}->{$gene}->{'pheno'}->{$pheno}->{'conf'}->{$conf}->{'no'} = $input;
    }
}

export_and_clear_stored_data();

sub export_and_clear_stored_data {
    # Don't bother with anything unless there are some basic required data: a variant name, a gene, and some phenotype.
    if (     ( exists $data_ref->{$curr_var_rnai}->{'header'} )
         and ( exists $data_ref->{$curr_var_rnai}->{'gene'}   )
         and ( exists $data_ref->{$curr_var_rnai}->{'pgene'}  ) ) {

        # Exactly *once*, print a top empty line.
        if (! $print_start ) { 
            print "\n";
            $print_start = 1;
        }

        # Print variation/RNAi title line:
        print "$data_ref->{$curr_var_rnai}->{'header'}\n";

        # One or more genes are affected by it:
        my @pre_g1_list = sort keys %{ $data_ref->{$curr_var_rnai}->{'gene'} };
        my @g1_list = grep { /WBGene\d+/ } @pre_g1_list;
        foreach my $g1 (@g1_list) { 

            # We only care about the ones with phenotypes:
            if ( exists $data_ref->{$curr_var_rnai}->{'pgene'}->{$g1} ) {

                # First, print each affected gene's name line:
                print "$data_ref->{$curr_var_rnai}->{'gene'}->{$g1}\n";

                # Then review the data for each possible phenotype.
                foreach my $p1 ( sort keys %{ $data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'} } ) {

                    # Require a source of evidence for each phenotype to avoid slop:
                    if ( exists $data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'}->{$p1}->{'conf'} ) { 

                        # For each evidence source...
                        foreach my $c1 ( sort keys %{ $data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'}->{$p1}->{'conf'} } ) { 

                            # Absolutely dichotomous choice.  Yes or no.  "No" overrides "yes!
                            if ( exists $data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'}->{$p1}->{'conf'}->{$c1}->{'no'} ) {
                                print "$data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'}->{$p1}->{'conf'}->{$c1}->{'no'}\n";
                            }
                            elsif ( exists $data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'}->{$p1}->{'conf'}->{$c1}->{'yes'} ) {
                                print "$data_ref->{$curr_var_rnai}->{'pgene'}->{$g1}->{'pheno'}->{$p1}->{'conf'}->{$c1}->{'yes'}\n";
                            }
                        }
                    }
                }
            }
        }
        # Print an empty line to separate each record.
        print "\n";
    }
    # Clean up the whole data structure, to keep RAM footprint light and avoid cross-talk between variants/RNAis:
    undef $data_ref->{$curr_var_rnai};
}

