#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

# Sample input:
#
# Position_Matrix : "Bartel_miRNA_5pFlank"
# Description    "Motif found upstream from microRNA genes."   Paper_evidence "WBPaper00024441" // pmid15317971
# Type           Frequency
# Background_model  A  0.314
# Background_model  C  0.186
# Background_model  G  0.186
# Background_model  T  0.314
# Site_values    A    0   1   0   0   0   0   0   0  16   4
# Site_values    C    31  8  33  33   4  31  33  23  14  22
# Site_values    G    0   0   0   0  27   2   0   0   0   0
# Site_values    T    2  24   0   0   2   0   0  10   3   7
# Sites_used   33
# Remark        "From paper WBPaper00024441/pmid15317971; consensus CTCCGCCC."

my $i         = 0;
my $matrix_id = q{};
my $brief_id  = q{};

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A Position_Matrix \s+ : \s+ \" ([^\"\s]+) \" \s* /xms ) { 
        $brief_id = $1;
        $brief_id = "Brief_id      \"$brief_id.pfm\"";

        $i++;
        $matrix_id = sprintf "%08i", $i;        

        my $header = 'Position_Matrix : "WBPmat' . $matrix_id . "\"";
        $data_ref->{'matrix'}->{$matrix_id}->{'header'}   = $header;
        $data_ref->{'matrix'}->{$matrix_id}->{'brief_id'} = $brief_id;
    } 
    elsif ( $input =~ /\A Description /xms ) { 
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
            if ( $input =~ /\A (Remark \s+ \" .+) [;] \s+ consensus \b .* \s (\S+) (\. \" .*) \z/xms ) { 
                my $remark_text1 = $1;
                my $consensus    = $2;
                my $remark_text2 = $3;

                $consensus = "Consensus     \"$consensus\"";
                $data_ref->{'matrix'}->{$matrix_id}->{'consensus'} = $consensus;
                my $remark = $remark_text1 . $remark_text2;
                push @{ $data_ref->{'matrix'}->{$matrix_id}->{'text_lines'} }, $remark;
            }
            else { 
                push @{ $data_ref->{'matrix'}->{$matrix_id}->{'text_lines'} }, $input;
            }
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

