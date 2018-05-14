#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my %type2abbrev = (
    Archaea     => 'Arch_',
    Arthropods  => 'Arth_',
    Bacteria    => 'Bact_',
    Eukaryotes  => 'Euka_',
    FIXME       => 'FIXME_',
    Metagenome  => 'Metg_',
    Metazoa     => 'Mzoa_',
    Nematodes   => 'Nema_',
    Vertebrates => 'Vert_',
    Viruses     => 'Viru_',
);

my @input_files = @ARGV;

my $fasta = pop @input_files;

foreach my $input_table (@input_files) {
    open my $TABLE, '<', $input_table;
    while (my $input = <$TABLE>) { 
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
            my $type    = $1;
            my $protein = $2;
            if (! exists $type2abbrev{$type}) {
                die "Can't abbreviate $type\n";
            }
            $data_ref->{'protein'}->{$protein}->{'type'}->{$type} = 1;
        }
        else {
            die "From input table $input_table, can't parse input line: $input\n";
        }
    }
    close $TABLE;
}

open my $FASTA, '<', $fasta;
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ /\A > (([^\s\/]+) \/ \d+ [-] \d+ \b .*) \z/xms ) { 
        my $header     = $1;
        my $protein    = $2;
        my $new_prefix = q{};
        my $output     = q{};
        if ( exists $data_ref->{'protein'}->{$protein}->{'type'} ) {
            my @types      = sort keys %{ $data_ref->{'protein'}->{$protein}->{'type'} };
            my $type_count = @types;

            if ( $type_count <= 0 ) {
                die "Failed to count types associated with: $input\n";
            }
            elsif ( $type_count == 1 ) {
                $new_prefix = $types[0];
            }
            elsif ( $type_count >= 2 ) {
                # None of these should coexist with any other type:
                if (    ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Archaea'}    )
                     or ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Bacteria'}   )
                     or ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Metagenome'} ) 
                     or ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Viruses'}    ) ) {
                    die "Inconsistent types (@types) associated with: $input\n";
                }
                elsif ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Arthropods'} ) { 
                    if (    ( exists $data_ref->{'protein'}->{ $protein }->{'type'}->{'Nematodes'}   ) 
                         or ( exists $data_ref->{'protein'}->{ $protein }->{'type'}->{'Vertebrates'} ) ) { 
                        die "Inconsistent types (@types) associated with: $input\n";
                    }
                    $new_prefix = 'Arthropods';
                }
                elsif ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Nematodes'} ) {
                    if (    ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Arthropods'}  )
                         or ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Vertebrates'} ) ) {
                        die "Inconsistent types (@types) associated with: $input\n";
                    }
                    $new_prefix = 'Nematodes';
                }
                elsif ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Vertebrates'} ) {
                    if (    ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Arthropods'} )
                         or ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Nematodes'}  ) ) {
                        die "Inconsistent types (@types) associated with: $input\n";
                    }
                    $new_prefix = 'Vertebrates';
                }
                elsif ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Metazoa'} ) {
                    $new_prefix = 'Metazoa';
                }
                elsif ( exists $data_ref->{'protein'}->{$protein}->{'type'}->{'Eukaryotes'} ) {
                    $new_prefix = 'Eukaryotes';
                }
                else {
                    # No matter how hard I try, it is impossible to get NCBI's TaxID dump to *perfectly* match the usage of UniProt.
                    # Thus, the best way to deal with these edge cases is to flag them, but warn about them loudly too, so that they
                    #    can be fixed one by one.
                    $new_prefix = 'FIXME';
                    warn "Somehow failed to parse type, so flagged with prefix 'FIXME_'; fix this by hand!: $input\n";
                } 
            }
        }
        else {
            $new_prefix = 'FIXME';
            warn "Failed to identify type of protein, so flagged with prefix 'FIXME_'; fix this by hand! in: $input\n";
        }
        $new_prefix = $type2abbrev{$new_prefix};
        $output     = '>' . $new_prefix . $header;
        print "$output\n";

    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Failed to parse FASTA header: $input\n";
    }
    else {
        print "$input\n";
    }
}

