#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

# Sample input:
# 
# Position_Matrix : "Tavazoie_motif_1"
# Description       "CTGCGTCT: motif number 1 predicted by Elemento and Tavazoie (2005)." Paper_evidence "WBPaper00033078"  // pmid15693947
# Type              Frequency
# Site_values       A 0 0 0 0 0 0 0 0
# Site_values       C 1 0 0 1 0 0 1 0
# Site_values       G 0 0 1 0 1 0 0 0
# Site_values       T 0 1 0 0 0 1 0 1
# Remark  "Evolutionarily conserved promoter motif, predicted in WBPaper00033078/pmid15693947; sense-strand orientation bias (5' to 3') with respect to promoter."

my $i         = 0;
my $matrix_id = q{};
my $brief_id  = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A Position_Matrix \s+ : \s+ \" ([^\"\s]+) \" \s* /xms ) { 
        $brief_id = $1;
        $brief_id =~ s/Tavazoie/Elemento_2005/;
        $brief_id = "Brief_id      \"$brief_id.pfm\"";

        $i++;
        $matrix_id = sprintf "%08i", $i;        

        my $header = 'Position_Matrix : "WBPmat' . $matrix_id . "\"";
        $data_ref->{'matrix'}->{$matrix_id}->{'header'}   = $header;
        $data_ref->{'matrix'}->{$matrix_id}->{'brief_id'} = $brief_id;
    } 
    elsif ( $input =~ /\A Description /xms ) { 
        if ( $input =~ /\A Description \s+ \" ([ACGT]+) [:] \s+ motif/xms ) { 
            my $consensus = $1;
            $consensus = "Consensus     \"$consensus\"";
            $data_ref->{'matrix'}->{$matrix_id}->{'consensus'} = $consensus;
            $input =~ s/Description\s+\"[ACGT]+[:]\s+motif/Description \"Motif/;
        }
        # Adjust spaces to make things line up more pleasingly in the reformatted version:
        $input =~ s/\ADescription\s+/Description   /;
        $data_ref->{'matrix'}->{$matrix_id}->{'description'} = $input;
    }
    elsif ( $input =~ /\A Type /xms ) { 
        if ( $input !~ /\A Type \s+ Frequency /xms ) { 
            die "Matrix $brief_id is not a PFM\n";
        }
        else { 
            # Again, adjust spaces to make things line up more pleasingly in the reformatted version:
            $input =~ s/\AType\s+/Type          /;
            $data_ref->{'matrix'}->{$matrix_id}->{'type'} = $input;
        }
    }
    else { 
        if ( $input =~ /\S/xms ) { 
            if ( $input =~ /\ASite_values/xms ) { 
                $input =~ s/\ASite_values\s+/Site_values   /;
            }
            if ( $input =~ /\ARemark/xms ) {
                $input =~ s/\ARemark\s+/Remark        /;
            }
            push @{ $data_ref->{'matrix'}->{$matrix_id}->{'text_lines'} }, $input;
        }
    }
}

my @matrices = sort keys %{ $data_ref->{'matrix'} };

print "\n";

foreach my $matrix_id1 (@matrices) { 
    print "$data_ref->{'matrix'}->{$matrix_id1}->{'header'}\n";

    print "$data_ref->{'matrix'}->{$matrix_id1}->{'brief_id'}\n";
    print "$data_ref->{'matrix'}->{$matrix_id1}->{'description'}\n";
    print "$data_ref->{'matrix'}->{$matrix_id1}->{'type'}\n";

    if ( exists $data_ref->{'matrix'}->{$matrix_id1}->{'consensus'} ) {
        print "$data_ref->{'matrix'}->{$matrix_id1}->{'consensus'}\n";
    }

    foreach my $matrix_text ( @{ $data_ref->{'matrix'}->{$matrix_id1}->{'text_lines'} } ) { 
        print "$matrix_text\n";
    }
    print "\n";
}

