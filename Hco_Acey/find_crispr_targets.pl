#!/usr/bin/env perl

# find_crispr_targets.pl -- Erich Schwarz <ems394@cornell.edu>, 5/2/2018.
# Purpose: read in 1+ FASTA DNA sequences; find possible CRISPR targets defined by variable criteria.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $data_ref;

my $dna       = q{};
my $header    = "DNA\tOri\tType\tStart\tSeq\tEnd\tPAM\tGCfrac\tInstances"; 
my $fasta     = ();
my $threshold = 'any';

# Keep the order of input sequences explicitly:
my @dnas = ();

my $help;

my %acceptable = (
    'all' => 1,
    'g1'  => 1,
    'gg2' => 1,
);

GetOptions ( 'fasta=s'     => \$fasta,
             'threshold=s' => \$threshold,
             'help'        => \$help, );

if ( (! $fasta) or (! $acceptable{$threshold} ) or $help ) { 
    die "Format: find_crispr_targets.pl\n",
        "    --fasta|-f       [FASTA DNA sequence file to search]\n",
        "    --threshold|-t   [stringency of site: 'gg2' (most stringent), 'g1' (less stringent), 'all' (default]\n",
        "    --help|-h        [print this message]\n",
        ;
}
    
open my $FASTA, '<', $fasta;
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) /xms ) {
        $dna = $1;
        if ( exists $data_ref->{'dna'}->{$dna} ) {
            die "Redundant input sequence: $dna\n";
        }
        else {
            push @dnas, $dna;
        }
    }
    elsif ( $input =~ /\S/xms ) {
        if (! $dna) {
            die "Sequence input undefined\n";
        }

        # delete blanks; reject non-DNA; make uppercase
        $input =~ s/\s//g;
        if ( $input =~ / [^ACGTacgt] /xms ) { 
            die "Sequence $dna has non-standard residues: $input\n";
        }
        $input =~ tr/acgt/ACGT/;
        $data_ref->{'dna'}->{$dna}->{'seq'} .= $input;
    }
}

foreach my $dna1 (@dnas) {
    my $seq = $data_ref->{'dna'}->{$dna1}->{'seq'};
    my @seq_nt = split //, $seq;

    my $len = length($seq);
    $len--;

    my $i = 0;
    my $j = 19;

    while ( $j <= $len) {
        my @site_nt = @seq_nt[$i..$j];

        my $site    = join q{}, @site_nt;
        my $revcomp = revcomp($site);

        my $cg_count = ( $site =~ tr/CG/CG/ );
        my $frac_cg  = ($cg_count / 20);

        my $start_nt = ($i + 1);
        my $end_nt   = ($j + 1);

        my $pam = q{};
        if ( ($j + 3) <= $len ) {
            my $k = ($j + 1);
            my $l = ($j + 3);
            my @poss_pam_nt  = @seq_nt[$k..$l];
            my $possible_pam = join q{}, @poss_pam_nt;
            if ( $possible_pam =~ /\A [ACGT] GG \z/xms ) {
                $pam = $possible_pam;
                $pam =~ tr/ACGT/acgt/;
            }
        }

        my $rev_pam = q{};
        if ( $i >= 3 ) {
            my $g = ($i - 3);
            my $h = ($i - 1);
            my @poss_rev_pam_nt  = @seq_nt[$g..$h];
            my $possible_rev_pam = join q{}, @poss_rev_pam_nt;
            $possible_rev_pam    = revcomp($possible_rev_pam);
            if ( $possible_rev_pam =~ /\A [ACGT] GG \z/xms ) { 
                $rev_pam = $possible_rev_pam;
                $rev_pam =~ tr/ACGT/acgt/;
            }
        }

        if ( ( $site =~ /[ACGT]{18} ([ACGT]{2})/xms ) and ($pam) ) {
            my $site_end_nt2 = $1;
            my $type         = q{};

            if ( $site_end_nt2 =~ /\A GG \z/xms ) {
                $type = 'gg2';
            }
            elsif ( $site_end_nt2 =~ /\A [ACGT] G \z/xms ) { 
                $type = 'g1'
            }
            else {
                $type = 'any';
            }

            print "$header\n" if $header;
            $header = q{};

            if (    ( $threshold eq 'any' )
                 or ( ( $threshold eq 'g1' )  and ( $type eq 'g1' ) )
                 or ( ( $threshold eq 'gg2' ) and ( $type eq 'gg2' ) )
               ) {

                # get the count of observed instances for this sequence:
                my $instances = 1;
                if ( exists $data_ref->{'target'}->{$site} ) {
                    $instances = $data_ref->{'target'}->{$site};
                }
                print "$dna1\t+\t$type\t$start_nt\t$site\t$end_nt\t$pam\t$frac_cg\t$instances\n";
                $instances++;
                $data_ref->{'target'}->{$site} = $instances;
            }
        }
        elsif ($pam) {
            die "Can't parse site: $dna1\t+\t[type?]\t$start_nt\t$site\t$end_nt\t$pam\t$frac_cg\n";
        }

        if ( ( $revcomp =~ /[ACGT]{18} ([ACGT]{2})/xms ) and ($rev_pam) ) {
            my $revcomp_end_nt2 = $1;
            my $type            = q{};

            if ( $revcomp_end_nt2 =~ /\A GG \z/xms ) {
                $type = 'gg2';
            }
            elsif ( $revcomp_end_nt2 =~ /\A [ACGT] G \z/xms ) {
                $type = 'g1'
            }
            else {
                $type = 'any';  
            }
            print "$header\n" if $header;
            $header = q{};
            if (    ( $threshold eq 'any' ) 
                 or ( ( $threshold eq 'g1' )  and ( ( $type eq 'g1' ) or ( $type eq 'gg2' ) ) ) 
                 or ( ( $threshold eq 'gg2' ) and ( $type eq 'gg2' ) )
               ) {

                # get the count of observed instances for this sequence:
                my $instances = 1;
                if ( exists $data_ref->{'target'}->{$revcomp} ) {  
                    $instances = $data_ref->{'target'}->{$revcomp};
                }
                print "$dna1\t-\t$type\t$end_nt\t$revcomp\t$start_nt\t$rev_pam\t$frac_cg\t$instances\n";
                $instances++;
                $data_ref->{'target'}->{$site} = $instances;
            }
        }
        elsif ($rev_pam) {
            die "Can't parse revcom site: $dna1\t-\t[type?]\t$end_nt\t$site\t$start_nt\t$rev_pam\t$frac_cg\n";
        }

        $i++;
        $j++;

    }
}
close $FASTA;

sub revcomp {
    my $in_string = $_[0];
    $in_string =~ tr/[acgtACGT]/[tgcaTGCA]/;
    my @in_residues = split //, $in_string;
    @in_residues = reverse @in_residues;
    my $out_string = join q{}, @in_residues;
    return $out_string;
}   
