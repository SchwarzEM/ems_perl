#!/usr/bin/env perl

# extract_genome_subset.pl -- Erich Schwarz <ems394@cornell.edu>, 10/17/2008; 
#     significant updates on 10/16/2011; major debugging by 10/22/2011!; generalized to any input sequences on 11/6/2013; fixed revcomp-on-[+] bug, 11/13/2013.
# Purpose: given DNA sequence FASTA and file of subset seq. coords., get subset FASTA from DNA., with the option of user-specified names or header info.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my %opts      = ();
my @chr_names = ();
my %chrs      = ();

GetOptions ( 'table=s' => \$opts{'table_file'},
             'fasta=s' => \$opts{'fasta_file'}, 
             'info'    => \$opts{'user_specified_header_info'},
             'names'   => \$opts{'user_specified_names'},
             'help'    => \$opts{'help'},        );

if (    $opts{'help'} 
     or (! $opts{'table_file'}                                              ) 
     or (! $opts{'fasta_file'}                                              ) 
     or ( $opts{'user_specified_header_info'} and $opts{'user_specified_names'} ) ) { 
    die "\n",
        "Format: extract_genome_subset.pl\n",
        "          --info|-i   [headers use block IDs, not names given in the input table (default)]\n",
        "          --names|-n  [headers use gene (etc.) names given in the coordinate table; mutually exclusive with --info|-i]\n",
        "          --table|-t  [table of coordinates of genomic blocks]\n",
        "                      [ chr_name  nt1  nt2  orient  gene(etc.)name ]\n",
        "                      [ tab-delimited; orient is + or - (i.e., FOR or REV); info becomes name or header text]\n",
        "          --fasta|-f  [genome FASTA to extract, perhaps already hard- or soft-masked]\n",
        "          --help      [print this message and quit]\n",
        "\n",
        ;
}

if ( (! exists $opts{'user_specified_header_info'} ) and (! exists $opts{'user_specified_names'} ) ) { 
    $opts{'user_specified_header_info'} = 1;
}

my $i    = q{};
my $chr  = q{};
my $nt1  = q{};
my $nt2  = q{};
my $info = q{};
my $ori  = q{};

my $data_ref;

# Note: this first gets the coordinates, then uploads *only*
# genomic residues that match those coordinates.  Trying things
# in the reverse order maxed-out my system's RAM...

open my $TABLE, '<', $opts{'table_file'} or die "Can't open interval table file $opts{'table_file'}: $!";
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\d+) \t (\d+) \t ([+]|[-]) \t ([^\t]+) \z /xms ) { 
        $chr  = $1;
        $nt1  = $2;
        $nt2  = $3;
        $ori  = $4;
        $info = $5;

        my $block_id = $chr . ':' . $nt1 . '..' . $nt2;
        # This shouldn't be necessary, but it's failsafe:
        if ( $nt1 > $nt2 ) { 
            ($nt1, $nt2) = ($nt2, $nt1);
            # Ignore $ori eq '+; -'.
            if ( $ori eq '+' ) { 
                $ori = '-';
            }
            elsif ( $ori eq '-' ) { 
                $ori = '+';
            }
        }

        # Store data for later use versus genome.
        # Note that $nt1 and $nt2 get recorded in their natural coordinates (1-based counting).
        # That is fine -- but when extracting sequences from a hashref using 0-based Perl-style 
        #     counting of residues, perform '--' before using these values to extract nucleotides!

        $data_ref->{'block_id'}->{$block_id}->{'chr'}       = $chr;
        $data_ref->{'block_id'}->{$block_id}->{'start_nt'}  = $nt1;
        $data_ref->{'block_id'}->{$block_id}->{'stop_nt'}   = $nt2;
        $data_ref->{'block_id'}->{$block_id}->{'ori'}       = $ori;
        $data_ref->{'block_id'}->{$block_id}->{'header'}    = $info;

        push @chr_names, $chr;
    }
    else { 
        die "From interval table file $opts{'table_file'}, cannot parse input line: $input\n";
    }
}
close $TABLE or die "Can't close filehandle to $opts{'table_file'}: $!";

# Do it this way to allow any chromosome names to be used, 
#     and to keep them in the same order as the input table.
@chr_names = uniq @chr_names;

open my $FASTA, '<', "$opts{'fasta_file'}" or die "Can't open FASTA file $opts{'fasta_file'}: $!";

while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ / \A > (\S+) /xms ) {
        $chr = $1;
        $i = -1;
        $chrs{$chr} = 1;
    }
    elsif ($input =~ /[a-zA-Z]/) {
        $input =~ s/[^a-zA-Z]//g;
        my @residues = split //, $input;
        foreach my $residue (@residues) { 
            $i++;
            $data_ref->{'genome_seq'}->{$chr}->{$i} = $residue;
        }
    }
}
close $FASTA or die "Can't close filehandle to $opts{'fasta_file'}: $!";

foreach my $block_id ( sort keys %{ $data_ref->{'block_id'} } ) { 
    $chr  = $data_ref->{'block_id'}->{$block_id}->{'chr'};
    $nt1  = $data_ref->{'block_id'}->{$block_id}->{'start_nt'};
    $nt2  = $data_ref->{'block_id'}->{$block_id}->{'stop_nt'};
    $ori  = $data_ref->{'block_id'}->{$block_id}->{'ori'};
    $info = $data_ref->{'block_id'}->{$block_id}->{'header'};

    # Enforce the reality of the demanded coordinates:
    if (! exists $chrs{$chr}) { 
         die "The tabular coordinates listed here do not, apparently, exist in the FASTA sequence:\n$chr\t$nt1\t$nt2\t$ori\t$info\n";
    }

    my $block_seq = q{};
    foreach my $nt ($nt1..$nt2) { 
        # Here, the number must shift from 1-based natural coordinates to 0-based Perl hashref coordinates:
        $nt--;
        if (! exists $data_ref->{'genome_seq'}->{$chr}->{$nt} ) {
            die "For block ID \"$block_id\" (header \"$info\", orientation \"$ori\") in chr. $chr", 
                " and nt range $nt1 to $nt2, failed to record sequence of chr. $chr, nt $nt!\n",
                ;
        }
        $block_seq .= $data_ref->{'genome_seq'}->{$chr}->{$nt};
    }
    # Revcomp only if specifically needed; ignore both '+' and '+; -'.
    if ( $data_ref->{'block_id'}->{$block_id}->{'ori'} eq '-' ) {
        $block_seq = revcomp($block_seq);
    }

    # Choice of previous hard-coded header format, or newer user-specified name format:
    if ( $opts{'user_specified_header_info'} ) {
        print ">$block_id  $info  $ori\n";
    }
    if ( $opts{'user_specified_names'} ) { 
        print ">$info  $ori  $block_id\n";
    }

    # FASTA format of seq. lines:
    my @output_lines 
        = unpack("a60" x (length($block_seq)/60 + 1), $block_seq);
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

