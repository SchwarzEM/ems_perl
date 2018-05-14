#!/usr/bin/env perl

# seqnumber_gffs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/11/2008.
# Purpose: add seqnumbers to GFF from FASTA.

use strict;
use warnings; 
use Getopt::Long;

my ( $gff_to_fix, $fasta_to_scan );

my $name      = q{};
my $count     = 0;
my %seq2size  = ();
my %seq_sized = ();

GetOptions ( 'gff:s'   => \$gff_to_fix, 
             'fasta:s' => \$fasta_to_scan, );

if ( (! $gff_to_fix) or (! $fasta_to_scan) or ($#ARGV > 1) ) { 
    die 'Format:',
        ' seqnumber_gffs.pl',
        ' --fasta/-f [FASTA] --gff/-g [GFF]',
        "\n",
        ;
}

open my $FASTA, '<', $fasta_to_scan 
    or die "Can't open FASTA $fasta_to_scan: $!";
while (my $input = <$FASTA>) { 
    chomp $input;
    if ($input =~ / \A > (\S+) /xms) { 
        my $new_name = $1;
        countseq($name,$count);
        $name = $new_name;
        $count = 0;
    }
    if (     ($name) 
         and ($input !~ / \A > /xms) 
         and ($input =~ / \w+  /xms) 
       ) { 
        $input =~ s/\s//g;
        # No error-check on FASTA input characters.
        $count += length($input);
    }
}
# Clear stored count at end-of-file.
countseq($name,$count);
close $FASTA or die "Can't close filehandle for $fasta_to_scan: $!";

open my $GFF, '<', $gff_to_fix 
    or die "Can't open GFF $gff_to_fix: $!";
while (my $input = <$GFF>) { 
    chomp $input;

    # Just print comment lines unchanged.
    if ($input =~ /\A \s* \# /xms) {
        print "$input\n";
    }

    # Process non-comment lines.
    if ($input !~ /\A \s* \# /xms) { 

        # Require sane input.
        if ( $input !~ / \A [^\t]+ \t (?: [^\t]* \t ){7} [^\t]* /xms ) { 
            die "Apparently misformatted: $input\n";
        }

        # When a GFF line shows up:
        if ( $input =~  /\A ([^\t]+) \t /xms ) { 
            $name = $1;

            # Prevent double-writing of Sequence lines:
            if ( $input =~ / \A 
                             [^\t]+   \t 
                             [^\t]+   \t
                             Sequence \t
                             1        \t
                             \d+      \t
                             (?: [^\t]+ \t ){3}
                             Sequence 
                             [ ] (?:\")? [^\"]+ (?:\")? 
                           /xms ) {
                die "$gff_to_fix already has:\n", "$input\n";
            }

            # Only try to add Seq. lines if FASTA was read:
            if ( (! exists $seq_sized{$name} ) and ( -e $fasta_to_scan ) ) { 

                # If, then, no seq. length available, die loudly:
                if (! exists $seq2size{$name}) {
                    my $swansong = "GFF $gff_to_fix" 
                                   . " has a sequence $name"
                                   . " that was not read"
                                   ;
                    if ($fasta_to_scan) { 
                        $swansong .= " from FASTA $fasta_to_scan";
                    }
                    die "$swansong\n";
                }

                # Otherwise, add a Sequence line, e.g.:
                # Chr1 . Sequence 1 14972282 . + . Sequence Chr1

                print "$name\t",
                      ".\t",
                      "Sequence\t",
                      "1\t",
                      "$seq2size{$name}\t",
                      ".\t",
                      "+\t",
                      ".\t",
                      "Sequence \"$name\"\n",
                      ;
                $seq_sized{$name} = 1;
            }
        }
        # Pass all other lines on normally.
        print "$input\n";
    }
}
close $GFF or die "Can't close filehandle for $gff_to_fix: $!";

sub countseq {
    ($name, $count) = @_;
    if ( ($count >= 1) and ($name =~ /\S/) ) { 
        $seq2size{$name} = $count;
    }
}

