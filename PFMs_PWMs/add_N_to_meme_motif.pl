#!/usr/bin/env perl

# add_N_to_meme_motif.pl -- Erich Schwarz <emsch@caltech.edu>, 11/10/2011.
# Purpose: given a MEME-style text motif with ACGT alphabet, create a version with ACGTN which can be used to scan masked genomic DNA.

use strict;
use warnings;

my $status = 1;
my $text1  = q{};
my $text2  = q{};
my $length = q{};
my $spacer = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A MEME [ ] version [ ] \d\S* .* \z /xms ) { 
        if ( $status != 1 ) { 
            die "Status is $status, so can't parse: $input\n";
        }
        else { 
            $status = 3;
            print "$input\n";
        }
    }
    elsif ( $input =~ /\A (ALPHABET= [ ] ACGT) (\s*) \z /xms ) { 
        $text1 = $1;
        $text2 = $2;
        if ( $status != 3 ) { 
            die "Status is $status, so can't parse: $input\n";
        }
        else { 
            $status = 4;
            print $text1, 'N', "$text2\n";
        }
    }
    elsif ( $input =~ /\A strands: [ ] \+ [ ] \- \s* \z/xms ) { 
        if ( $status != 4 ) {
            die "Status is $status, so can't parse: $input\n";
        }
        else {
            $status = 5;
            print "$input\n";
        }
    }
    elsif ( $input =~ /\A Background [ ] letter [ ] frequencies [ ] \( from [ ] .+ \) : \s* \z/xms ) {
        if ( $status != 5 ) {
            die "Status is $status, so can't parse: $input\n";
        }
        else {
            $status = 6;
            print "$input\n";
        }
    }
    elsif ( $input =~ /\A \s* (A \s \d\.(\d+) (\s) C \s \d\.\d+ \s  G \s \d\.\d+ \s T \s \d\.\d+) \s* \z/xms ) {
        $text1  = $1;
        $length = $2;
        $spacer = $3;
        $length = length($length);
        my $N_val = '0' x $length;
        $N_val = '0.' . $N_val;

        if ( $status != 6 ) {
            die "Status is $status, so can't parse: $input\n";   
        }
        else {
            $status = 7;
            print $text1, $spacer, 'N', $spacer, "$N_val\n";
        }
    }
    elsif ( $input =~ /\A MOTIF [ ]+ \S .* \z/xms ) { 
        # For multi-motif files, this line can happen 
        #    either after a full set of information or after the very first line,
        #    so allow two different values of $status here:
        #    the first value is single-motif/file, the second only happens with recurrent motifs.
        if ( ( $status != 7 ) and ( $status != 2 ) ) {
            die "Status is $status, so can't parse: $input\n";  
        }
        else {
            $status = 8;
            print "$input\n";
        }
    }
    elsif ( $input =~ /\A BL \s+ MOTIF \s+ \S+ \s+ width=0 \s+ seqs=0 \s* \z/xms ) { 
        if ( $status != 8 ) {
            die "Status is $status, so can't parse: $input\n";
        }
        else {
            # Leave status unchanged, since this line appears to be entirely optional.
            print "$input\n";
        }
    }
    elsif ( $input =~ /\A (letter-probability [ ] matrix: [ ] alength= [ ]) 4 ([ ] w= [ ] \d+ [ ] nsites= [ ] \d+ [ ] E= [ ] \S+ \s*) \z/xms ) {
        $text1 = $1;
        $text2 = $2;
        if ( $status != 8 ) {
            die "Status is $status, so can't parse: $input\n";
        }
        else {
            $status = 9;
            print $text1, '5', "$text2\n";
        }
    }
    elsif ( $input =~ /\A ( \s* \d\.(\d+) (\s+) \d\.\d+ \s+ \d\.\d+ \s+ \d\.\d+ ) \s* \z/xms ) {
        $text1 = $1;
        $length = $2;
        $spacer = $3;
        $length = length($length);
        my $N_val = '0' x $length;
        $N_val = '0.' . $N_val;

        if ( $status != 9 ) {
            die "Status is $status, so can't parse: $input\n";
        }
        else {
            # Leave status unchanged, since this sort of line will appear many times.
            print $text1, $spacer, "$N_val\n";
        }
    }
    elsif ( $input =~ /\A \s* \z/xms ) { 
        if ( $status == 9 ) { 
            $status = 2;
        }
        print "$input\n";
    }
    else { 
        die "Status is $status; can't parse: $input\n";
    }
}

