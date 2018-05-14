#!/usr/bin/env perl

use strict;
use warnings;

my $orth_table = $ARGV[0];
my $data_table = $ARGV[1];

my %names2merged = ();

open my $ORTHS, '<', $orth_table or die "Can't open ortholog table $orth_table: $!";
while (my $input = <$ORTHS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $name1 = $1;
        my $name2 = $2;
        if ( $name1 eq $name2 ) {
            if ( $name1 eq 'Gene' ) {
                $names2merged{'Gene'} = 'Gene';
            }
            else { 
                die "Unacceptable line from ortholog table $orth_table: $input\n";
            }
        } 
        else { 
            my $merged = $name1 . q{|} . $name2;
            $names2merged{$name1} = $merged;
            $names2merged{$name2} = $merged;
        }
    }
    else { 
        die "From ortholog table $orth_table, can't parse input: $input\n";
    }
}
close $ORTHS or die "Can't close filehandle to ortholog table $orth_table: $!";

open my $DATA, '<', $data_table or die "Can't open data table $data_table: $!";
while (my $input = <$DATA>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) (\t .*) \z/xms ) { 
        my $name = $1;
        my $text = $2;
        if ( exists $names2merged{$name} ) {
            my $output = $names2merged{$name} . $text;
            print "$output\n";
        }
    }
    else { 
        die "From data table $data_table, can't parse: $input\n";
    }
}
close $DATA or die "Can't close filehandle to data table $data_table: $!";

