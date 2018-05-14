#!/usr/bin/env perl

use strict;
use warnings;

my $pdb_names   = $ARGV[0];
my $broken_text = $ARGV[1];

my %tag2name = ();

open my $PDB, '<', $pdb_names or die "Can't open PDB names file $pdb_names: $!";
while (my $input = <$PDB>) {
    chomp $input;
    # E.g.: >gi|116666565|pdb|1VS5|R
    if ( $input =~ /\A > ( gi \| \d+ \| pdb \| ([A-Z0-9]+ \| [A-Z0-9]*) ) /xms ) {
        my $full_name = $1;
        my $tag       = $2;
        $tag =~ s/\|//;
        if ( exists $tag2name{$tag} ) { 
            die "Failure! There is an ambiguous tag, $tag, associated with the sequence name $full_name!\n";
        }
        $tag2name{$tag} = $full_name;
    }
    else { 
        die "Can't parse text from PDB names file $pdb_names: $input\n";
    }
}
close $PDB or die "Can't close filehandle to PDB names file $pdb_names: $!";

open my $BROKEN, '<', $broken_text or die "Can't open broken text file $broken_text: $!";
while (my $input = <$BROKEN>) {
    chomp $input;
    if ( $input =~ /\ASubsidiary [ ] database \t Accession [ ] no. \t Motif [ ] name \t PDB [ ] sequence\z/xms  ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A ([^\t]+ \t [^\t]+ \t [^\t]+) \t ([^\t]+) \z/xms ) {
        my $text1 = $1;
        my $tag   = $2; 
        if (! exists $tag2name{$tag} ) {
            die "Failure! There is an ambiguous tag, $tag, which cannot be parsed in broken text file $broken_text: $input\n";
        }
        my $full_name = $tag2name{$tag};
        print "$text1\t$full_name\n";
    } 
    else { 
        die "Failure!  Cannot correct text in broken text file $broken_text: $input\n";
    }
}
close $BROKEN or die "Can't close filehandle to broken text file $broken_text: $!";

