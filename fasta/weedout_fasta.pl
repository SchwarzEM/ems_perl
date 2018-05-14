#!/usr/bin/perl

# weeded_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/5/2007.
# Purpose: derivative of uniform_fasta.pl that censors zero-length, no M1, or int-stopped proteins.

use strict;
use warnings;

my $output_line     = q{};
my @output_lines    = ();
my $seq_name        = q{};
my %seqs2headers    = ();
my %sequences       = ();
my %errors_list_ref = ();

# Record the incoming FASTA data.

while (my $input_line = <>) { 
    chomp $input_line;
    if ($input_line =~ /\A > ( (\S+) .*) /xms) { 
        $seq_name = $2; 
        $seqs2headers{$seq_name} = $1;
        $sequences{$seq_name} = q{};
    }
    elsif ( $input_line =~ /\S/xms ) { 
        $sequences{$seq_name} .= $input_line;
    }
}

# Weed out bad sequences; print out the good ones.

foreach my $seq_name2 (sort keys %sequences) { 
    # Record as errors: no start 'M' codon; empty sequence (!); internal stop '*' codon(s).
    if ( $sequences{$seq_name2} =~ /[a-zA-Z]\*[a-zA-Z]/xms ) { 
        $errors_list_ref{$seq_name2}->{'internal stop codon(s)'} = 1; 
    }
    if ( $sequences{$seq_name2} !~ / \A (m|M) /xms ) { 
        $errors_list_ref{$seq_name2}->{'no starting M codon'} = 1;
    }
    if ( $sequences{$seq_name2} !~ / [a-zA-Z] /xms ) {
        $errors_list_ref{$seq_name2}->{'zero-length sequence'} = 1;
    }

    # Print all sequence without above bloopers.
    if ( ( $sequences{$seq_name2} =~ /[a-zA-Z]/xms ) 
       and ( $sequences{$seq_name2} !~ /[a-zA-Z]\*[a-zA-Z]/xms ) 
       and ( $sequences{$seq_name2} =~ / \A (m|M) /xms ) ) { 
        print ">$seqs2headers{$seq_name2}\n";
        @output_lines 
            = unpack("a60" x (length($sequences{$seq_name2})/60 + 1), $sequences{$seq_name2});
        foreach $output_line (@output_lines) { 
            if ($output_line =~ /\S/) { 
                print "$output_line\n";
            }
        }
    }
}

# Record erroneous sequences from the input FASTA.

my $date = join('.', &get_local_date());
my $errors = 'errors' . '.' . $date;
my $i = 0;
while (-e $errors) { 
    $i++;
    $errors =~ s/\-\d+\z//;
    $errors = $errors . "-$i";
    if ($i == 100_000) { 
        die "Dude!  Why are there 100,000 look-alike errors files?\n";
    }
}

open my $ERRORS, '>', $errors or die "Can't open error log $errors: $!";

foreach my $seq_name3 (sort keys %errors_list_ref) { 
    my @errors = sort keys %{ $errors_list_ref{$seq_name3} };
    my $error_line = join '; ', @errors;
    print {$ERRORS} "$seq_name3: $error_line\n";
}

close $ERRORS or die "Can't close filehandle for $errors: $!";

sub get_local_date { 
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;
}

