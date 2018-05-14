#!/usr/bin/env perl

# make_revcomp_seqs.pl -- Erich Schwarz <ems394@cornell.edu>, 11/13/2013.
# Purpose: given a DNA FASTA file, revcomp it while keeping the order of sequence names intact -- optionally with no indication of [REVCOMP] status (useful for fixing files that were mistakenly revcomp in the first place!).

use strict;
use warnings;

use strict;
use warnings;
use Getopt::Long;

my $fasta  = q{};
my $silent = q{};
my $seq_id = q{};
my $help;

my $data_ref;

GetOptions ( 'fasta=s' => \$fasta, 
             'silent'  => \$silent,
             'help'    => \$help,  );

if ( $help or (! $fasta) ) { 
        die "Format: make_revcomp_seqs.pl\n",
            "    --fasta|-f   [DNA FASTA sequence to revcomp, perhaps hard- or soft-masked]\n",
            "    --silent|-s  [do not append \" [REVCOMP]\" to header lines (the default action)]\n",
            "    --help       [print this message and quit]\n",
            ;
}

open my $FASTA, '<', "$fasta" or die "Can't open FASTA file $fasta: $!";
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ / \A > ((\S+) .*) \z/xms ) {
        my $header = $1;
        $seq_id    = $2;
        if ( exists $data_ref->{'seq'}->{$seq_id} ) { 
            die "In FASTA file $fasta, redundant sequence name: $seq_id\n";
        }
        $data_ref->{'seq'}->{$seq_id}->{'header'} = $header;
        push @{ $data_ref->{'seqlist'} }, $seq_id;
    }
    elsif ($input =~ /[a-zA-Z]/) {
        $input =~ s/[^a-zA-Z]//g;
        $data_ref->{'seq'}->{$seq_id}->{'sequence'} .= $input;
    }
}
close $FASTA or die "Can't close filehandle to $fasta: $!";

foreach my $seq_id1 (@{ $data_ref->{'seqlist'} }) { 
    my $header = $data_ref->{'seq'}->{$seq_id1}->{'header'};
    if (! $silent) { 
        $header .= ' [REVCOMP]';
    }
    print ">$header\n";

    my $sequence = $data_ref->{'seq'}->{$seq_id1}->{'sequence'};
    $sequence = revcomp($sequence);

    # FASTA format of seq. lines:
    my @output_lines 
        = unpack("a60" x (length($sequence)/60 + 1), $sequence);
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
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

