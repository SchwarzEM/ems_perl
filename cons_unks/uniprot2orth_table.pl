#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw{uniq};

my $data_ref;

my @accs = ();

my %ortho_class = (
    eggNOG => 1,
    GeneTree => 1,
    InParanoid => 1,
    OMA => 1,
    OrthoDB => 1,
    PhylomeDB => 1,
);

while ( my $input = <>) {
    chomp $input;

    # For each ID, zero out the ACCs array.
    if ( $input =~ /A ID \s+ \S+ /xms ) {
        @accs = ();
    }

    # For each line detecting ACCs associated with an ID, add them to whatever list exists.
    # There should normally be just one line, but that can't be guaranteed, given UniProt irregular formatting.
    elsif ( $input =~ /\A AC \s+ \S+ \s* /xms ) {

        # A single UniProt ID can have two or more UniProt accession numbers.
        if ( $input =~ /\A AC \s+ (\S+ ; .*) \z/xms ) { 
            my $acc_txt = $1;
            $acc_txt =~ s/;\s*\z//;
            my @new_accs = split /; /, $acc_txt;
            push @accs, @new_accs;
        }
        else { 
            die "Cannot parse input: $input\n";
        }
    }

    elsif ( $input =~ /\A DR \s+ (\S+) ; \s+ (\S+) ; /xms ) {
        my $db     = $1;
        my $db_acc = $2;
        if ( ( exists $ortho_class{$db} ) and ( @accs ) ) { 
            @accs = uniq @accs;  # This shouldn't be necessary, but is a failsafe against redundancy.
            foreach my $acc (@accs) {
                print "$db|$db_acc\t$acc\n";
            }
        }
    }
}

