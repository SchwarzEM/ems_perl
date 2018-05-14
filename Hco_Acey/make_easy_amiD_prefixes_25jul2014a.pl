#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %prot2type = ();

my %type2abbrev = (
    Archaea => 'Arch_',
    Arthropods => 'Arth_',
    Metagenome => 'Metg_',
    Metazoa => 'Mzoa_',
    Nematodes => 'Nema_',
    Vertebrates => 'Vert_',
    Viruses => 'Vir_',
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
            $prot2type{$protein} = $type2abbrev{$type};
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
    if ( $input =~ /\A > ((([^\s\/]+) \/ \d+) [-] \d+ \b .*) \z/xms ) { 
        my $header     = $1;
        my $seqname    = $2;
        my $protein    = $3;
        my $new_prefix = 'Bact_';  # default
        if (exists $prot2type{$protein}) { 
            $new_prefix = $prot2type{$protein};
        }
        my $output = '>' . $new_prefix . $header;
        print "$output\n";
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse FASTA header: $input\n";
    }
    else {
        print "$input\n";
    }
}

