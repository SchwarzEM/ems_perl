#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $reject_text   = q{};
my @reject_names  = ();
my %not_OK       = ();

my $header        = q{};
my $drugs_text    = q{};
my @drug_names    = ();
my @ok_drug_names = ();

my $print_fasta  = 0;

my %opts = ();

GetOptions ( 'fasta=s'  => \$opts{'input_fasta'}, 
             'reject=s' => \$opts{'reject_fasta'},
             'help'     => \$opts{'help'},         );

if ( $opts{'help'} or (! $opts{'input_fasta'} ) or (! $opts{'reject_fasta'} ) ) { 
    die "Format: prune_drugbank_targets.pl --fasta|-f [FASTA to filter] --reject|-r [FASTA with rejected drugs] --help|-h\n";
}

# Sample inputs:
# >drugbank_target|15 Voltage-dependent T-type calcium channel subunit alpha-1I (DB01388)
# >drugbank_target|20 Prostaglandin G/H synthase 1 (DB03783; DB04817)

open my $REJECT, '<', $opts{'reject_fasta'} or die "Can't open reject-list FASTA file $opts{'reject_fasta'}: $!";
while (my $input = <$REJECT>) { 
    chomp $input;
    if ( $input =~ / \A 
                     > 
                     .+ 
                     \( 
                     (.+) 
                     \) 
                     \s* 
                     \z /xms ) {
        $reject_text = $1;
        @reject_names = split /; /, $reject_text;
        foreach my $rname (@reject_names) { 
            $not_OK{$rname} = 1;
        }
    }
    elsif ( $input =~ / \A > /xms ) { 
        die "Can't parse input FASTA header from reject-list FASTA file $opts{'reject_fasta'}:\n $input\n";
    }
}
close $REJECT or die "Can't close filehandle to reject-list FASTA file $opts{'reject_fasta'}: $!";

open my $FASTA, '<', $opts{'input_fasta'} or die "Can't open FASTA file to be filtered, $opts{'input_fasta'}: $!";
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( ( $input !~ / \> /xms ) and $print_fasta ) { 
        print "$input\n";
    }
    elsif ( $input =~ / \> /xms ) {
        $print_fasta = 0;
        if ( $input =~ / \A
                     (>.+)
                     \(
                     (.+)
                     \)
                     \s*
                     \z /xms ) {
            $header        = $1;
            $drugs_text    = $2;
            @drug_names    = split /; /, $drugs_text;
            # Re-initialize this:
            @ok_drug_names = ();
            foreach my $dname (sort @drug_names) { 
                if (! exists $not_OK{$dname} ) { 
                    push @ok_drug_names, $dname;
                }
            }
        }
        if (@ok_drug_names) { 
            $print_fasta = 1;
            my $drugs_text = join '; ', @ok_drug_names;
            print "$header($drugs_text)\n";
        }
    }
    elsif ( $input =~ / \A > /xms ) {
        die "Can't parse input FASTA header from filtered FASTA file $opts{'input_fasta'}:\n $input\n";
    }
}
close $FASTA or die "Can't close filehandle to FASTA file to be filtered, $opts{'input_fasta'}: $!";


