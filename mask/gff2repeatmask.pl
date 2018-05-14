#!/usr/bin/env perl

# gff2repeatmask.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/20/2008.
# Purpose: hard-/soft-mask complex/simple repeats; req. FASTA and GFF2.

# Caveats:
#
# 1. This script presupposes correct FASTA format.
# 2. It will obliterate previous softmasking!

use strict;
use warnings;
use Getopt::Long;

my $fasta_file;
my $gff_file;
my $mask_type;
my $repeats;

my $i              = q{};
my $seq_name       = q{};
my $start_nt       = q{};
my $stop_nt        = q{};
my %mask_ranges    = ();

GetOptions ( 'gff:s'       => \$gff_file,
             'fasta:s'     => \$fasta_file,
             'mask_type:s' => \$mask_type, 
             'repeats:s'   => \$repeats,    );

$gff_file   ||= q{};
$fasta_file ||= q{};
$mask_type  ||= 'soft';
$repeats    ||= 'both';

if (       (! -e $gff_file                               ) 
        or (! -e $fasta_file                             )
        or ( $mask_type !~ /\A(soft|hard)\z/             ) 
        or ( $repeats !~ /\A(both|complex|simple)\z/ ) ) { 
    die 'Format: gff2repeatmask.pl ',
        ' --gff/-g=[GFF2 of repeats]',
        ' --fasta/-f=[FASTA]',
        ' --mask_type/-m=[soft*|hard]',
        ' --repeats/-r=[both*|complex|simple]  (*default)',
        "\n",
    ;
}

open my $GFF, "<", "$gff_file"
    or die "Can't open GFF file $gff_file: $!";

my $pattern1 = '(\S+) \s+ tandem       \s+'
               . ' tandem_repeat   \s+ (\d+) \s+ (\d+) \t';
my $pattern2 = '(\S+) \s+ inverted     \s+'
               . ' inverted_repeat \s+ (\d+) \s+ (\d+) \t';
my $pattern3 = '(\S+) \s+ RepeatMasker \s+'
               . ' repeat_region   \s+ (\d+) \s+ (\d+) \t';

while (my $input = <$GFF>) {
    chomp $input;

# sample acceptable input (tab-delimited) lines:
# I       tandem  tandem_repeat   18142   18171   53      .       .       Note "2 copies of 15mer"
# I       inverted        inverted_repeat 21551   21635   92      .       .       Note "loop 470, 5 gaps"
# I       RepeatMasker    repeat_region   845     1061    1823    .       .       Target "Motif:Ce000122" 68 284

    # Default, mask all repeats; optionally, only complex or simple ones:
    if  (    ( ( $repeats ne 'complex') 
                 and ( $input =~ / \A $pattern1 /xms ) )
          or ( ( $repeats ne 'complex') 
                 and ( $input =~ / \A $pattern2 /xms ) )
          or ( ( $repeats ne 'simple') 
                 and ( $input =~ / \A $pattern3 /xms ) ) ) {
        $seq_name = $1;

        # Convert from 1+ DNA counting to 0+ Perl-array counting:
        $start_nt = ($2 - 1);  

        # Ditto:
        $stop_nt  = ($3 - 1);

        # Do NOT use simple '##' in regex, which clobbers whole script!
        if ( $seq_name !~ / \A \#\# /xms ) { 
            foreach my $nt ($start_nt..$stop_nt) { 
                $mask_ranges{$seq_name}->{$nt} = 1;
            }
        }
    }
}
close $GFF;

open my $FASTA, '<', "$fasta_file"
    or die "Can't open FASTA file $fasta_file: $!";

# Should cope with huge FASTAs, by handling only one line at a time:
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ / \A > (\S+) /xms ) {
        $seq_name = $1;
        # First residue is '0':
        $i = 0;
        if (! exists $mask_ranges{$seq_name} ) {
            warn "No masking data were provided for $seq_name!\n";
        }
        print "$input\n";
    }
    elsif ($input =~ /[a-zA-Z]/) {
        # Reject incorrect residues:
        $input =~ s/[^a-zA-Z]//g;

        # Force uppercase before softmasking:
        $input =~ tr/[a-z]/[A-Z]/;

        my $output = q{};
        my @residues = split //, $input;
        foreach my $residue (@residues) { 
            # Only alter masked residues:
            if ( exists $mask_ranges{$seq_name}->{$i} ) { 
                if ( $mask_type eq 'soft' ) {
                    $residue =~ tr/[A-Z]/[a-z]/;
                }
                if ( $mask_type eq 'hard' ) {
                    $residue = 'N';
                }
            }

            # Do this even if nothing was altered:
            $output .= $residue;

            # Increment $i after reading residue '0', etc.:
            $i++;
        }
        # Lines come out with unchanged lengths + res. coords.:
        print "$output\n";
    }
}
close $FASTA
    or die "Can't close filehandle to $fasta_file: $!";

