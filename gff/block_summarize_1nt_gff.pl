#!/usr/bin/env perl

# block_summarize_1nt_gff.pl -- Erich Schwarz <ems394@cornell.edu>, 6/4/2018.
# Purpose: given a genome sequence and a 1-nt/line GFF3 annotation (e.g., 1-nt liftover), identify blocks of annotated and non-annotated nt.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $dna = q{};
my $gff = q{};
my $pos = 'annot';
my $neg = 'NON_ANN';

my @seqs    = ();
my $seq     = q{};
my $seq_len = 0;

my $data_ref;

my $help;

GetOptions ( 'dna=s' => \$dna,
             'gff=s' => \$gff,
             'pos=s' => \$pos,
             'neg=s' => \$neg,
             'help'  => \$help, );

if ( $help 
     or (! -e $dna     ) 
     or (! -e $gff     ) 
     or ( $pos eq $neg ) 
   ) { 
    die "Format: block_summarize_1nt_gff.pl\n",
        "    --dna|-d   [input DNA sequence being annotated]\n",
        "    --gff|-g   [1-nt/line GFF3 annotation]\n",
        "    --pos|-p   [positive annotation for positive nt; default \"annot\"]\n",
        "    --neg|-n   [negative annotation for negative nt; default \"NON_ANN\"]\n",
        "    [--pos and --neg may not be identical]\n",
        "    --help|-h  [print this message]\n",
        ;
}

open my $DNA, '<', $dna;
while (my $input = <$DNA>) {
    chomp $input;

    if ( $input =~ /\A > (\S+) /xms ) {
        my $new_seq = $1;

	if ( exists $data_ref->{'seq'}->{$new_seq} ) {
            die "Redundant sequence: $seq\n";
        }

        if ($seq) {
            if ( $seq_len < 1 ) {
                die "Sequence $seq has length of $seq_len\n";
            }
            $data_ref->{'seq'}->{$seq}->{'length'} = $seq_len;
            foreach my $i (1..$seq_len) {
                $data_ref->{'seq'}->{$seq}->{'nt'}->{$i}->{'status'} = $neg;
            }
            push @seqs, $seq;
        }

        $seq     = $new_seq;
        $seq_len = 0;
    }

    elsif ( $input =~ /\S/xms ) {
        if (! $seq ) {
            die "No sequence available for: $input\n";
        }

        $input =~ s/\s//g;
        if ( $input =~ /[^ACGTNacgtn]/xms ) { 
            die "Non-nucleotide residues for sequence $seq at: $input\n";
        }
        my $nt_len = length($input);
        $seq_len = ($seq_len + $nt_len);
    }
}

# Record remaining stored seq. data after done reading FASTA.
$data_ref->{'seq'}->{$seq}->{'length'} = $seq_len;
foreach my $i (1..$seq_len) {
    $data_ref->{'seq'}->{$seq}->{'nt'}->{$i}->{'status'} = $neg;
}
push @seqs, $seq;

# And zero out:
$seq     = q{};
$seq_len = 0;

close $DNA;

open my $GFF, '<', $gff;
while (my $input = <$GFF>) {
    chomp $input;

    # instance input line:
    # chrI_pilon      c_elegans.PRJNA13758.WS264.genomic.fa   nt      1835    1835    .       .       .       seq=G;orig_nt=I:1

    if ( $input =~ /\A (\S+) \t (\S+) \t nt \t (\d+) \t (\d+) \t \. \t \. \t \. \t seq=\S; orig_nt=(\S+:\d+) \z/xms ) {
        $seq         = $1;
        my $source   = $2;
        my $nt1      = $3;
        my $nt2      = $4;
        my $prev_nt  = $5;

        if ( $nt1 != $nt2 ) {
            die "Unequal nt positions in $input\n";
        }

        $data_ref->{'seq'}->{$seq}->{'nt'}->{$nt1}->{'status'}  = $pos;
        $data_ref->{'seq'}->{$seq}->{'nt'}->{$nt1}->{'source'}  = $source;
        $data_ref->{'seq'}->{$seq}->{'nt'}->{$nt1}->{'prev_nt'} = $prev_nt;
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
close $GFF;
$seq = q{};

foreach	my $seq1 (@seqs) {
    $seq_len = $data_ref->{'seq'}->{$seq1}->{'length'};

    my $block_status = q{};
    my $block_start  = 0;
    my $block_end    = 0;
    my $new_status   = 0;

    foreach my $nt (1..$seq_len) {
        $new_status = $data_ref->{'seq'}->{$seq1}->{'nt'}->{$nt}->{'status'};

        # Automatically begins a new block when previous $block_status is null.
        if ( $new_status ne $block_status) { 
            if ($block_status) {
                if ( $block_start < 1 ) {
                    die "For sequence $seq1, bad $block_status block start: $block_start\n";
                }
                if ( $block_end < $block_start ) {
                    die "For sequence $seq1, bad $block_status block end: $block_end\n";
                }
                $data_ref->{'seq'}->{$seq1}->{'block_start'}->{$block_start}->{'block_status'} = $block_status;
                $data_ref->{'seq'}->{$seq1}->{'block_start'}->{$block_start}->{'block_end'} = $block_end;
            }
            $block_status = $new_status;
            $block_start  = $nt;
            $block_end    = $nt;
        }
        else {
            # Let *end* +1 nt with each new line; leave status and start fixed.
            $block_end    = $nt;
        }
    }
    # Record final stored data at end of sequence.
    $data_ref->{'seq'}->{$seq1}->{'block_start'}->{$block_start}->{'block_status'} = $block_status;
    $data_ref->{'seq'}->{$seq1}->{'block_start'}->{$block_start}->{'block_end'} = $block_end;
}

foreach my $seq2 (@seqs) {
    my @block_starts = sort { $a <=> $b } keys %{ $data_ref->{'seq'}->{$seq2}->{'block_start'} };
    foreach my $block_start (@block_starts) {
        my $block_status = $data_ref->{'seq'}->{$seq2}->{'block_start'}->{$block_start}->{'block_status'};
        my $block_end    = $data_ref->{'seq'}->{$seq2}->{'block_start'}->{$block_start}->{'block_end'};

        my $prev_block_start = $data_ref->{'seq'}->{$seq}->{'nt'}->{$block_start}->{'prev_nt'};
        my $prev_block_end   = $data_ref->{'seq'}->{$seq}->{'nt'}->{$block_end}->{'prev_nt'};

        my $source = $data_ref->{'seq'}->{$seq}->{'nt'}->{$block_start}->{'source'};
        if ( $source ne $data_ref->{'seq'}->{$seq}->{'nt'}->{$block_end}->{'source'} ) {
            warn "Discordant sources for block in $seq2 (nt $block_start vs. nt $block_end)\n";
        }

        print "$seq2\t";
        print "$source\t";
        print "$block_status\t";
        print "$block_start\t";
        print "$block_end\t";
        print ".\t.\t.\n";
        print 'orig_nt', $prev_block_start, '-', $prev_block_end;
        print "\n";
    }
}

