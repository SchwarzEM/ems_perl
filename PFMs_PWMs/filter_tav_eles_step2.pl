#!/usr/bin/env perl

# filter_tav_eles_step2.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/31/2010.
# Purpose: Given a kludge subset of Tavazoie predictions, make a rough but usable .ace-like output.

use strict;
use warnings;
use TFBS::PatternGen::SimplePFM;

my $serial_no = q{};
my $nt_seq    = q{};
my $ori1      = q{};
my $ori2      = q{};
my $go_term   = q{};
my $go_pval   = q{};  
my $transfac  = q{};

my $blank_hdr = 1;

my @sequences = ();
my @remarks   = ();

# Read in sequence lines from raw list or FASTA:
while (my $input = <>) { 
    chomp $input;

# Sample input lines:
# 
#     *... first interval is *not* simple \t!
# 2        CTGCGTCTC      381.9   557.5   0       2000    3.8e-09**       1.0e+00 3.81                                                    
# 3        AGACGCAGA      320.0   630     0       2000    1.0e+00 2.2e-08**       3.41                                                    
# 4        CGACACTCC      241.5   234     0       1500    4.2e-01 5.8e-01 2.49            positive regulation of growth   8.95e-08        

    if ( $input =~ / \A (\d+) 
                        \s+ 
                        ([ACGT]+) 
                        \t (?:[^\t]*\t){4}  
                        (\S+) 
                        \t 
                        (\S+) 
                        [^\t]* \t [^\t]* \t [^\t]* \t
                        ([^\t]*)
                        \t
                        ([^\t]*)
                        \t
                        ([^\t]*)
                        /xms ) { 
        $serial_no = $1;
        $nt_seq    = $2; 
        $ori1      = $3;
        $ori2      = $4;
        $go_term   = $5;
        $go_pval   = $6;
        $transfac  = $7;

        if ( $ori2 =~ /\*\*\z/xms ) { 
            $nt_seq = revcomp($nt_seq);
            ($ori1, $ori2) = ($ori2, $ori1);
        }
        # Specifically reverse GATA back, so that it remains visibly GATA -- as done in the paper:
        if ( $serial_no == 5 ) { 
            $nt_seq = revcomp($nt_seq);
            ($ori1, $ori2) = ($ori2, $ori1);
        }
        push @sequences, $nt_seq;

        # First time only, print a filler space:
        if ($blank_hdr) { 
            print "\n";
            $blank_hdr = 0;
        }

        # For each line, print a .ace text:
        print q{Position_Matrix : "Tavazoie_motif_}, 
              $serial_no, 
              q{"}, 
              "\n", 
              ;
        print "Description       \"$nt_seq: motif number $serial_no predicted by Elemento and Tavazoie (2005).\"",
              " Paper_evidence \"WBPaper00033078\"  // pmid15693947\n",
              ;
        print "Type              Frequency\n";
    }


    # Generate, format, and print frequency lines:
    my $pfm = TFBS::PatternGen::SimplePFM->new(-seq_list=>\@sequences);
    my @countlines = split /\n/, ($pfm->pattern->prettyprint);

    foreach my $line (@countlines) {
        $line =~ s/(\[|\])//g;
        $line =~ s/\A[ ]+//;
        $line =~ s/[ ]+\z//;
        $line =~ s/[ ]{2,}/ /g;
        print "Site_values       $line\n";
    }
 
    if ( ( $ori1 =~ /\*\*\z/xms )  or ( $ori2 =~ /\*\*\z/xms ) or $go_term or $go_pval or $transfac ) { 
        @remarks = ();
        push @remarks, "Evolutionarily conserved promoter motif, predicted in WBPaper00033078/pmid15693947";
        if ( $ori1 =~ /\*\*\z/xms ) { 
            push @remarks, "sense-strand orientation bias (5' to 3') with respect to promoter;";
        }
        if ( $ori2 =~ /\*\*\z/xms ) {
            push @remarks, "antisense-strand orientation bias (3' to 5') with respect to promoter;";
        }
        if ( ( $go_term =~ /\w/xms ) and ( $go_pval =~ /\w/xms ) ) {
            push @remarks, "GO term: $go_term (p-value, $go_pval)";
        }
        if ($transfac) { 
            push @remarks, "Transfac matches: $transfac";
        }
        @remarks = grep { $_ =~ /\w/xms } @remarks;
        my $remark_line = join q{; }, @remarks;
        $remark_line =~ s/(;){2,}/;/g;
        $remark_line =~ s/\t/ /g;
        $remark_line =~ s/\s+,/,/g;
        $remark_line =~ s/\A\s+//;
        $remark_line =~ s/[ ]+/ /g;
        $remark_line =~ s/;\z//;
        if ( $remark_line =~ /\w/xms ) { 
            print "Remark  \"";
            print "$remark_line.\"\n";
        }
    }
    print "\n";
    $serial_no = q{};
    $nt_seq    = q{};
    $ori1      = q{};
    $ori2      = q{};
    $go_term   = q{};
    $go_pval   = q{};
    $transfac  = q{};
    @sequences = ();
    @remarks   = ();
}

sub revcomp { 
    my $input_sequence = $_[0];
    $input_sequence = reverse($input_sequence);
    if ( $input_sequence =~ /[^acgtACGT]/xms ) { 
        die "Not currently designed to parse any",
            " letters but a, c, g, t, A, C, G, or T!\n",
            ;
    }
    $input_sequence =~ tr/acgtACGT/tgcaTGCA/;
    return $input_sequence;
}

