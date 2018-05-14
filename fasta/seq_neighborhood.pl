#!/usr/bin/env perl

# seq_neighborhood.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/10/2010.
# Purpose: get neighborhoods of 1+ nt around defined words in a FASTA sequence; optionally, look for revcomp as well (good for DNA!), and/or trim flank 'N/n' residues.

use strict;
use warnings;
use Getopt::Long;

my $center     = q{};
my $revcenter  = q{};
my $flank_size = 0;
my $revcomp;
my $trim_Nn_or_Xx;
my $help;

my $seq_name     = q{};
my $header       = q{};
my @infiles      = ();
my @output_lines = ();
my $output_line  = q{};
my $i            = 0;
my $data_ref;

GetOptions ( 'fasta=s{,}' => \@infiles,
             'center=s'   => \$center,
             'size=i'     => \$flank_size,
             'revcomp'    => \$revcomp,
             'trim'       => \$trim_Nn_or_Xx,
             'help'       => \$help,     );

if ( $help or (! @infiles) or (! $center) or (! $flank_size ) ) {
    die "Format: seq_neighborhood.pl\n",
        "            --fasta|-f   [input FASTA file(s)]\n",
        "            --center|-c  [sequence at center of neighborhood]\n",
        "            --size|-s    [non-zero flank size in residues]\n",
        "            --revcomp|-r [optionally, check for both center seq. and its\n",
        "                             reverse-complement (useful for DNA, useless for protein!)\n",
        "            --trim|-t    [optionally, trim neighborhoods so that they keep no N/n or X/x residues]\n",
        "            --help|-h    [print this message]\n",
        ;
}

foreach my $infile1 (@infiles) { 
    open my $INFILE, '<', $infile1 or die "Can't open input FASTA file $infile1: $!";
    while (my $input = <$INFILE>) { 
        chomp $input;
        if ($input =~ /\A > ( (\S+) .*) /xms) { 
            $header   = $1;
            $seq_name = $2;
            $data_ref->{'input_seq'}->{$seq_name}->{'header'} = $header;
        }
        elsif ( $input =~ /\A > /xms ) { 
            die "Aberrant input line: $input\n";
        }
        elsif ( $input =~ /[a-zA-Z]/xms ) { 
            $data_ref->{'input_seq'}->{$seq_name}->{'sequence'} .= $input;
        }
    }
    close $INFILE or die "Can't close filehandle to input FASTA file $infile1: $!";
}

foreach my $seq_name2 ( sort keys %{ $data_ref->{'input_seq'} } ) {
    $i = 0;
    scan_seq_for_wordzones($seq_name2, $center, $flank_size);
    if ($revcomp) {
        $revcenter = revcomp($center);
        scan_seq_for_wordzones($seq_name2, $revcenter, $flank_size);
    }
}

foreach my $subseq ( sort keys %{ $data_ref->{'output_seq'} } ) { 
    print ">$subseq    $data_ref->{'output_seq'}->{$subseq}->{'header'}\n";
    @output_lines
        = unpack("a60" x (length($data_ref->{'output_seq'}->{$subseq}->{'sequence'})/60 + 1), $data_ref->{'output_seq'}->{$subseq}->{'sequence'});
    foreach $output_line (@output_lines) {
        if ($output_line =~ /\S/) {
            print "$output_line\n";
        }
    }
}

sub scan_seq_for_wordzones {
    my $_parent_seq  = $_[0];
    my $_center      = $_[1];
    my $_flank_size  = $_[2];   
    my $_subseq      = q{};
    my $_subseq_name = q{};
    while ( $data_ref->{'input_seq'}->{$_parent_seq}->{'sequence'} 
                =~ / ( [A-Za-z]{$_flank_size} $_center [A-Za-z]{$_flank_size} ) /xmsg ) {
            $_subseq = $1;
            $i++;
            $_subseq_name = $_parent_seq . '_' . $i;
            if ( exists $data_ref->{'output_seq'}->{$_subseq_name} ) {
                die "Redundant subsequence name from $_parent_seq: $_subseq_name\n";
            }
            if ($trim_Nn_or_Xx) { 
                if ( $_subseq =~ / \A .*? [NnXx]+ ([^NnXx]* $_center .*) \z/xms ) { 
                    $_subseq = $1;
                }
                if ( $_subseq =~ / \A (.* $_center [^NnXx]*) [NnXx]+ .*? \z/xms ) { 
                    $_subseq = $1;
                }
            }
            $data_ref->{'output_seq'}->{$_subseq_name}->{'sequence'} = $_subseq;
            $data_ref->{'output_seq'}->{$_subseq_name}->{'header'} 
                = $data_ref->{'input_seq'}->{$_parent_seq}->{'header'};
    }    
}

sub revcomp {
    my $in_string = $_[0];
    $in_string =~ tr/[acgtACGT]/[tgcaTGCA]/;
    my @in_residues = split //, $in_string;
    @in_residues = reverse @in_residues;
    my $out_string = join q{}, @in_residues;
    return $out_string;
}   

