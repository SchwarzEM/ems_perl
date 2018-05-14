#!/usr/bin/env perl

# DIANA-1.1.2.pl

# originally by Melissa Michelitsch,
# with modifications by Erich Schwarz <emsch@its.caltech.edu>, 1/18/2002 and 8/5/2008.
#
# Derived from the original DIANA, by:
#   Melissa Michelitsch <mmichel@itsa.ucsf.edu>
#       [e-mail address on 1/18/02]
#   Melissa Michelitsch <melissa_michelitsch@chiron.com>
#       [e-mail address after UCSF address stops being used]
# 
# Version 1.0 of DIANA was described in the article:
# 
# Michelitsch and Weissman (2000).  A census of 
# glutamine/asparagine-rich regions: implications for
# their conserved function and the prediction of novel prions.
# Proc. Natl. Acad. Sci. U.S.A. vol. 97, pp. 11910-11915.
# 
# DIANA-1.1.pl:   1/18/2002.  Changes to make DIANA easier to run.
# DIANA-1.1.1.pl: 1/30/2002.  Bugfix, to correctly parse FASTA files.
# DIANA-1.1.2.pl: 1/30/2002.  Change to make hitlists simpler/more useful.
# DIANA-2.0.1.pl: 9/5/2008.   Major, major code revision from 2002 to 2008 practices.

use strict;
use warnings;

# accept a protein FASTA file as flow input.

my $protein_fasta_input = q{};
if ( -e $ARGV[0] ) { 
    $protein_fasta_input = $ARGV[0];
}
if ( ( $protein_fasta_input !~ /\A\w/xms ) 
     or ( $protein_fasta_input !~ /\w\z/xms ) ) { 
    $protein_fasta_input = 'protein_input';
}

my $seq_name               = q{};
my $diana_output           = $protein_fasta_input . ".diana-output";
my $filtered_diana_output  = $protein_fasta_input . ".diana-filtered-output";
my %diana_seq              = ();
my $diana_hitlist          = $protein_fasta_input . ".diana-hitlist";

open my $DIANA_OUTPUT, '>', $diana_output
    or die "Can't open DIANA output file $diana_output. $! \n";
open $FILTERED_DIANA_OUTPUT, '>', $filtered_diana_output"
    or die "Can't open filtered DIANA output file $filtered_diana_output! $! \n";
open $DIANA_HITLIST, '>', $diana_hitlist
    or die "Can't open DIANA hitlist file $diana_hitlist! $! \n";

open my $PROTEIN_FASTA_INPUT, '<', $protein_fasta_input
    or die "Can't open input file $protein_fasta_input. $!\n";

while (my $input = <$PROTEIN_FASTA_INPUT>) {
    chomp $input;
    if (    ( $input =~ / \A > gi\| (\w+) \| /xms )
         or ( $input =~ / \A > (\S+) \s+     /xms )     
         or ( $input =~ / \A > (\S+) \z      /xms ) ) { 
        $seq_name = $1;
    }
    elsif ( $input =~ / \A > /xms ) { 
        die "Malformatted input line: $input\n";
    }
    else { 
        $input =~ s/\s//g;
        $diana_seq{$seq_name} .= $input;
    }
}
close $PROTEIN_FASTA_INPUT 
    or die "Can't close filehandle to $protein_fasta_input. $!\n";

my $seqnum = 1; 

my $count_of_seqs = keys %diana_seq;

# my @array = split /\*/, $foobar; 

##Make Letter-sequence array ## 

foreach my $seq_name2 (sort keys %diana_seq) { 
    my $x = 0;
    my $y = 0;
    my $z = 0;
    my $r = 0;
    my @letarray = split //, $diana_seq{$seq_name2};
    my $ltrcount = length $diana_seq{$seq_name2};

    my $maxQ = 0;

    if ($ltrcount <= 79) {
        $x = $ltrcount;
    } 
    else {
        $x = 79;
    } 

##make 80-mers## 

    $ATmer = join q{}, @letarray[$z..$x]; 

    @mer = split /Q|N|q|n/, $ATmer;

    if ($letarray[$x] !~ /Q|N|q|n/) {
        $numQ = @mer - 1; 
    }
    else {
        $numQ = @mer; 
    }

## ref. seq. by # of Q's ## 

    while ($letarray[$x] ne q{}) {
        if ($z != 0) { 
            if ( $letarray[$x] =~ /Q|N|q|n/ ) { 
                $numQ++;
            } 
            $k = ($z - 1);

            if ( $letarray[$k] =~ /Q|N|q|n/ ) {
                $numQ--;
            }
        }

##determine max Q ## 

        if ( $maxQ <= $numQ ) { 
            $ATmer = join q{}, @letarray[$z..$x]; 
            $Qarray{$numQ} = $ATmer; 
            $maxQ = $numQ; 
            ##write to assoc. array by max number Q calling up assoc. seq. ## 
        }
        $x++; 
        $z++;
        $r++;
    }

    $seqstr = $Qarray{$maxQ}; 
    @Qcnt   = split /Q|q/, $seqstr; 
    $numQ   = (@Qcnt - 1); 
    @Ncnt   = split /N|n/, $seqstr;
    $numN   = (@Ncnt - 1);

##Output##

    my $output_line = $title[$b] 
                      . "\t" 
                      . $Qarray{$maxQ} 
                      . "\t" 
                      . $ltrcount 
                      . "\t" 
                      . $maxQ 
                      . "\t" 
                      . $numQ 
                      . "\t" 
                      . $numN 
                      ;
    print $DIANA_OUTPUT "$output_line\n";
    if ( $maxQ > 29 ) { 
        print $FILTERED_DIANA_OUTPUT "$output_line\n";
        print $DIANA_HITLIST "$title[$b]\n";
    }
} 

close $DIANA_OUTPUT;
close $FILTERED_DIANA_OUTPUT;
close $DIANA_HITLIST;

