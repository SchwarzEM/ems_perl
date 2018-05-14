#!/usr/bin/perl

# report_step1.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/19/2007.
# Purpose: step 1 of churning out Tables 1-2 for 2007 Hox paper.

use strict;
use warnings;

# Formatting constants for sprintf in output:
my $X = 0;
my $Y = 0;

my $fasta   = q{};
my %fa_lens = ();

while (my $input = <>) { 
    chomp $input;
    if ($input =~ / \A > (\S+) /xms) { 
        $fasta = $1;
        if (length($fasta) > $X) { 
            $X = length($fasta);
        }
        if ($fa_lens{$fasta}) { 
            die "Multiple entries for $fasta\n";
        }
        $fa_lens{$fasta} = 0;
    }
    elsif ( ( $input =~ /\S/ ) and ($fasta) ) { 
        $input =~ s/\s//g;
        $fa_lens{$fasta} += length($input);
        if ( length($fa_lens{$fasta}) > $Y) { 
            $Y = length($fa_lens{$fasta});
        }
    }
}

# Add some whitespace ...?
$Y += $Y + ($Y % 3);  # len. of added commas, apparently all I need!

foreach my $seq (sort keys %fa_lens) { 

    # Number of residues (nt, aa) in $seq:
    my $seqlen = q{};

    # Print extra line before starting new fosmid's data:
    if ($seq =~ /\. tfa \z/xms ) {
        print "\n";
    }

    # Don't change $seq in place! make print-friendly copy.
    my $seq_to_print = sprintf "%".'-'.$X."s", $seq;

    # Fixed-width left-justified field, ergo no tab.
    print "$seq_to_print";

    # OMIT tab-delimiting, for now.
    # print "\t";

    # Weird regex-magic from Perl Cookbook 2.16, p. 84:
    $seqlen = commify($fa_lens{$seq});

    # Append residue type based solely on suffix -- kludge!
    if ($seq =~ /\. tfa \z/xms ) {
        $seqlen .= " nt";
    }
    if ($seq !~ /\. tfa \z/xms ) {
        $seqlen .=  " aa";
    }

    # Fixed-width format and print:
    $seqlen = sprintf "%".$Y."s", $seqlen;
    print $seqlen;

    # End the line:
    print "\n";
}

# Perl Cookbook 2.16, p. 84 -- baffling, but it works:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

