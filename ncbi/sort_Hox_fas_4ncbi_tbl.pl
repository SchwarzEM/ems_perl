#!/usr/bin/env perl

# sort_HoxPaper_fastas.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/29/2008.
# Purpose: sort-out into:
#     one free-standing FASTA for each DNA contig
#     a single multi-FASTA monolith for all CDSes, complete or not, from a given contig

use strict;
use warnings;

if ($#ARGV != 1) { 
    die "Format: ./sort_HoxPaper_fastas.pl [contigs file] [proteins file]\n";
}

my $contig2prots_ref;

my $contig_file  = $ARGV[0];
my $protein_file = $ARGV[1];

# Put here for ease of code revision; note the escape key for regex!
my $ncbi_prefix  = 'lcl\|';

open my $CONTIG_FILE, '<', $contig_file 
   or die "Can't open contig file $contig_file: $!";

my $contig_sequence = q{};

while (my $input = <$CONTIG_FILE>) { 
    if ( ( $input !~ / \A > $ncbi_prefix \S+ \s .* \n \z/xms )
         and ( $input =~ / \A > /xms ) ) { 
        die "Can't parse DNA contig FASTA header: $input";
    }
    if ( $input =~ /\A \> $ncbi_prefix (\S+) \s .* \n \z/xms ) { 
        $contig_sequence = $1;
    }
    if ( ( $contig_sequence =~ /\S/xms ) and ( $input =~ /\S/xms ) ) { 
        $contig2prots_ref->{$contig_sequence}->{'dna'} .= $input;
    }
}

close $CONTIG_FILE or die "Cannot close filehandle to $contig_file: $!";
        
open my $PROTEIN_FILE, '<', $protein_file
    or die "Can't open contig file $protein_file: $!";

my $protein_sequence = q{};
my $parent_contig    = q{};

while (my $input = <$PROTEIN_FILE>) {
    if ( ( $input !~ /\A > $ncbi_prefix \S+ \s .* \n \z/xms ) 
         and ( $input =~ /\A >/xms ) ) { 
        die "Can't parse protein FASTA header: $input";
    }
    if ( $input =~ /\A \> $ncbi_prefix (\S+) \s .* \[gene=$ncbi_prefix(\S+)\] .* \n \z/xms ) {
        $parent_contig    = $1;
        $protein_sequence = $2;
        if (! exists $contig2prots_ref->{$parent_contig}->{'dna'} ) { 
            die "Whoa, cowboy! protein $protein_sequence",
                " has no recognizable parent contig.\n",
                ;
        }
    }
    if ( ( $protein_sequence =~ /\S/xms ) and ( $input =~ /\S/xms ) ) {            
        $contig2prots_ref->{$parent_contig}->{'protein'} .= $input;
    }
}    

close $PROTEIN_FILE or die "Cannot close filehandle to $protein_file: $!";

foreach my $contig (sort keys %{ $contig2prots_ref } ) { 
    open my $OUTPUT, '>', "$contig.fsa"
        or die "Cannot open DNA output file $contig.fsa: $!";
    print {$OUTPUT} $contig2prots_ref->{$contig}->{'dna'};
    close $OUTPUT 
        or die "Cannot close filehandle to DNA output file $contig: $!";
    
    open $OUTPUT, '>', "$contig.pep"
        or die "Cannot open full-protein output file $contig.pep: $!";
    print {$OUTPUT} $contig2prots_ref->{$contig}->{'protein'};
    close $OUTPUT 
        or die "Cannot close filehandle to full-protein output file $contig.prot.fa: $!";
}

