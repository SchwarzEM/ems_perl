#!/usr/bin/env perl

# sort_fasta.pl -- Erich Schwarz <emsch@caltech.edu>, 11/9/2011.
# Purpose: given one or more FASTA files and an ordered list of either sequence names or partial sequence names (e.g., genes), create a FASTA file which follows the list order.

use strict;
use warnings;
use Getopt::Long;

my @fastas = ();
my $FASTA;
my $list   = q{};

my $header       = q{};
my $name         = q{};
my $i            = 0;
my @output_lines = ();

my $seqnames;
my $partnames;
my $data_ref;
my $help;

GetOptions ( "fasta=s{,}" => \@fastas,
             "list=s"     => \$list,
             "seqnames"   => \$seqnames,
             "partnames"  => \$partnames,
             "help"       => \$help, );

if ( $help 
     or (! @fastas                         ) 
     or (! $list                           ) 
     or ( (! $seqnames) and (! $partnames) ) 
     or ( $seqnames and $partnames         ) ) {
    die "sort_fasta.pl\n",
        "    --fasta|-f [1+ FASTA file(s)]      --list|-l [list of either seq. or part-seq. (gene) names]\n",
        "    --seqnames|-s <or> --partnames|-p  [list has either exact names, or subnames such as gene IDs, but not both]\n", 
        "    --help|-h\n",
        ;
}

foreach my $fasta (@fastas) { 
    open $FASTA, '<', $fasta or die "Can't open FASTA input file $fasta: $!";
    while (my $input = <$FASTA>) { 
        chomp $input;
        if ( $input =~ /\A > ((\S+) .*) \z/xms ) { 
            $header = $1;
            $name   = $2;
            $i++;
            if ( exists $data_ref->{'seqname'}->{$name} ) { 
                die "Redundant sequence, $name, in input FASTA(s): @fastas\n";
            }
            $data_ref->{'seqname'}->{$name}->{'header'}     = $header;
            $data_ref->{'seqname'}->{$name}->{'order_seen'} = $i;
        }
        elsif ( $input =~ /\A > /xms ) { 
            die "From FASTA file $fasta, can't parse input line: $input\n";
        }
        elsif ( $input =~ / \S+ /xms ) { 
            $input =~ s/\s//g;
            $data_ref->{'seqname'}->{$name}->{'sequence'} .= $input;
        }
    }
    close $FASTA or die "Can't close filehandle to FASTA input file $fasta: $!";
}

# Reset $i for the serialization of names in $list.
$i = 0;

open my $LIST, '<', $list or die "Can't open list of names $list: $!";
while (my $input = <$LIST>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \z /xms ) { 
        $name = $1;
        $i++;
        $data_ref->{'listname'}->{$name}->{'order_listed'} = $i;
    }
    else { 
        die "Did not parse list input line: \"$input\"\n";
    }
}
close $LIST or die "Can't close filehandle to list of names $list: $!";

foreach my $seqname (grep { /\S/ } sort keys %{ $data_ref->{'seqname'} } ) { 
    if (! exists $data_ref->{'seqname'}->{$seqname}->{'order_seen'} ) { 
        die "Failed to get initial order for sequence $seqname!\n";
    }
    foreach my $listname (sort keys %{ $data_ref->{'listname'} } ) { 
        if (! exists $data_ref->{'listname'}->{$listname}->{'order_listed'} ) { 
            die "Failed to get order for listed name $listname!\n";
        }

        # Simple case where we have an exact order list for each exact sequence name.
        elsif ($seqnames) { 
            if (! exists $data_ref->{'listname'}->{$seqname}->{'order_listed'} ) {
                die "Failed to get list order for sequence name $seqname!\n";
            }
            $data_ref->{'seqname'}->{$seqname}->{'final_listing'} 
                = $data_ref->{'listname'}->{$seqname}->{'order_listed'};
        }

        # More interesting situation where we have partial names -- e.g., WBGene IDs -- that are ordered.
        elsif ($partnames) { 
            # In some cases, /xms could be tricky for getting a subseq, but it shouldn't be here.
            if ( $seqname =~ /$listname/xms ) { 
                if (   (exists $data_ref->{'seqname'}->{$seqname}->{'final_listing'}    ) 
                   and ( $data_ref->{'seqname'}->{$seqname}->{'final_listing'} 
                             != $data_ref->{'listname'}->{$listname}->{'order_listed'} ) ) { 
                    die "Inconsistency: the list name order gives two different",
                        " rankings for sequence $seqname,",
                        " $data_ref->{'seqname'}->{$seqname}->{'order_listed'}",
                        " versus  $data_ref->{'listname'}->{$listname}->{'order_listed'}!\n",
                        ;
                }
                $data_ref->{'seqname'}->{$seqname}->{'final_listing'} 
                    = $data_ref->{'listname'}->{$listname}->{'order_listed'};
            }
        }
    }
}

# Get the final list of sequences which I will both order and export.
my @final_seq_list = ();

# Only do this with genes which got a 'final_listing' value.
#    I.e., that were wanted by being included in the list!

foreach my $candidate_seq ( keys %{ $data_ref->{'seqname'} } ) { 
    if ( exists $data_ref->{'seqname'}->{$candidate_seq}->{'final_listing'} ) { 
        push @final_seq_list, $candidate_seq;
    }
}

# Have the order seen as a fallback baseline ordering,
#    with list order being imposed on top of it.
#    This helps keep a sensible order in situations where
#    a partial name list (e.g. genes) is used to order
#    multiple (1+/gene) sequences associated with genes.

@final_seq_list = sort { $data_ref->{'seqname'}->{$a}->{'final_listing'} <=> $data_ref->{'seqname'}->{$b}->{'final_listing'} } 
                  sort { $data_ref->{'seqname'}->{$a}->{'order_seen'}    <=> $data_ref->{'seqname'}->{$b}->{'order_seen'}   }
                  @final_seq_list ;

foreach my $final_seq (@final_seq_list) { 
    print ">$data_ref->{'seqname'}->{$final_seq}->{'header'}\n";
    @output_lines 
        = unpack("a60" x (length($data_ref->{'seqname'}->{$final_seq}->{'sequence'})/60 + 1), 
                 $data_ref->{'seqname'}->{$final_seq}->{'sequence'});
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

