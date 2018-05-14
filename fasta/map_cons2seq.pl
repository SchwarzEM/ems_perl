#!/usr/bin/perl

# map_id2seq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/11/2008.
# Purpose: w/ 2 aligned FASTAs, mark identities of index to consensus seq. as upper-case of index.

use strict;
use warnings;
use Getopt::Long;

my $index_file     = q{};
my @index_residues = ();

my $cons_file      = q{};
my @cons_residues  = ();

my %seq2header     = ();

my $index_CONSeq   = q{};

GetOptions ( "index=s"     => \$index_file,
             "consensus=s" => \$cons_file, ); 

seq2array_w_header($index_file, \@index_residues);
seq2array_w_header($cons_file, \@cons_residues);

sub seq2array_w_header { 
    my ($infile, $array_ref) = @_;
    open my $FH, '<', $infile 
         or die "Can't open sequence file $infile: $!";
    my @residues       = ();
    while (my $input = <$FH>) { 
        chomp $input;
        if ($input =~ /\A > (\S+ .*) \z /xms) { 
            my $header = $1;
            $seq2header{$infile} = $header;
        }
        if ( ($input !~ /\A > \S+ /xms)
           and ($input =~ /[a-zA-Z\-]/xms) ) {
            $input =~ s/[^a-zA-Z\-]//g;
            my @residues = split //, $input;
            push @{ $array_ref }, @residues;
        }
    }
    close $FH or die "Can't close filehandle to $infile: $!";
}

foreach my $i (0..$#index_residues) { 
    my $residue = $index_residues[$i];
    if ( $residue ne '-' ) { 
        if ( $residue eq $cons_residues[$i] ) { 
            $residue =~ tr/[a-z]/[A-Z]/;
        }
        else { 
            $residue =~ tr/[A-Z]/[a-z]/;
        }
        $index_CONSeq .= $residue;
    }
}

print '>', $seq2header{$index_file}, "\n";
# print "$index_CONSeq\n";

my @output_lines 
   = unpack("a60" x (length($index_CONSeq)/60 + 1), $index_CONSeq);
foreach my $output_line (@output_lines) { 
    if ($output_line =~ /\S/) { 
        print "$output_line\n";
    }
}

