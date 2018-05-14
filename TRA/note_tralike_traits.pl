#!/usr/bin/env perl

# tralike_traits.pl -- Erich Schwarz <emsh@its.caltech.edu>, 10/6/2009.
# Purpose: given a FASTA file of protein sequences and args, extract TRA-like traits.

use strict;
use warnings;
use Getopt::Long;

my $seqid       = q{};
my ($lowest_len, $highest_len, $middle_len, $cterm_len);
my $seq_data_ref;

GetOptions ( 'lowest_len:i'  => \$lowest_len,
             'highest_len:i' => \$highest_len,
             'middle_len:i'  => \$middle_len,
             'cterm_len:i'   => \$cterm_len,  );

# Defaults for any variables not entered as arguments:

$lowest_len  ||= 160;    # Default: was 158 == (197 aa of D. melanogaster TRA * 0.8).
$highest_len ||= 515;    # Default: was 515 == (429 aa of C. capitata TRA * 1.2).
$middle_len  ||= 100;    # Really, the N-terminal flank of the C-terminus.
$cterm_len   ||= 60;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $new_seqid = $1;
        if ($seqid) { 
            process_rawseq( $seqid,       $lowest_len, 
                            $highest_len, $middle_len, 
                            $cterm_len,   $seq_data_ref, ); 
        }
        $seqid = $new_seqid;
    }
    elsif ( $input =~ / \S /xms ) { 
        $input =~ s/\s//g;
        $seq_data_ref->{$seqid}->{'rawseq'} .= $input;
    }
}
if ($seqid) {
    process_rawseq( $seqid,       $lowest_len, 
                    $highest_len, $middle_len, 
                    $cterm_len,   $seq_data_ref, ); 
}

print "Seq.\tLen.\tP_no\tP_rat\tRS_no\tR+S_no\n";

foreach my $seq_name ( sort keys %{ $seq_data_ref } ) { 
    print "$seq_name";
    print "\t";
    print "$seq_data_ref->{$seq_name}->{'length'}"; 
    print "\t";
    print "$seq_data_ref->{$seq_name}->{'P_count'}";
    print "\t";
    print "$seq_data_ref->{$seq_name}->{'P_ratio'}";
    print "\t";
    print "$seq_data_ref->{$seq_name}->{'RS_count'}";
    print "\t";
    print "$seq_data_ref->{$seq_name}->{'R_and_S_count'}";
    print "\n";
}

sub process_rawseq { 
    my $_seq_name     = $_[0];
    my $_lowest_len   = $_[1];
    my $_highest_len  = $_[2];
    my $_middle_len   = $_[3];
    my $_cterm_len    = $_[4];
    my $_data_ref     = $_[5];

    my $_raw_seq = $_data_ref->{$_seq_name}->{'rawseq'};
    my $_length = length($_raw_seq);
    if ( $_length < 1 ) { 
        die "Illegal data value: seq. $_seq_name; len. $_length!\n";
    }
    if ( ( $_length >= $_lowest_len ) and ( $_length <= $_highest_len ) ) { 
        # For $_start1 and $_start2: no '+ 1' because of 0-based residue counting.
        my $_start1        = $_length - $_cterm_len - $_middle_len;  
        my $_start2        = $_length - $_cterm_len;  # Solely for clarity;
        my $_non_cterm_len = $_length - $_cterm_len;  # these are the same!

        my $_non_cterm_seq = substr($_raw_seq, 0,       $_non_cterm_len);
        my $_mid_seq       = substr($_raw_seq, $_start1, $_middle_len);
        my $_cterm_seq     = substr($_raw_seq, $_start2, $_cterm_len);

        $_data_ref->{$_seq_name}->{'length'}        = $_length;
        my $_C_term_Pcount                           = ($_cterm_seq     =~ s/P/P/ig);
        my $_nonC_term_Pcount                        = ($_non_cterm_seq =~ s/P/P/ig);
        $_data_ref->{$_seq_name}->{'P_count'}       = $_C_term_Pcount;
        $_data_ref->{$_seq_name}->{'P_ratio'}       = ( $_C_term_Pcount / ($_nonC_term_Pcount + 0.001) );
        $_data_ref->{$_seq_name}->{'RS_count'}      = ($_mid_seq =~ s/RS/RS/ig);
        $_data_ref->{$_seq_name}->{'R_and_S_count'} = ( ($_mid_seq =~ s/R/R/ig) + ($_mid_seq =~ s/S/S/ig) );
    }
    if ( ( $_length < $_lowest_len ) or ( $_length > $_highest_len ) ) { 
        delete $_data_ref->{$_seq_name};  # Keep load on RAM low; make readout easy later.
    }
    return;
}

