#!/usr/bin/env perl

# wb_vs_yang_contigs_04sep2009.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/2/2009.
# Former version of: wb_vs_yang_contigs.pl
# Purpose: compare specified nt of WormBase (WS204) contigs to Yang's contig frags.

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my $WB;
my $wb = q{};

my $YANG_DIR;
my $yang_dir = 0;

my $YANG_FILE;
my $yang_file = q{};

my $prefix    = q{};

GetOptions ( "wb=s"     => \$wb,
             "yang=s"   => \$yang_dir,
             "prefix=s" => \$prefix, );

if ( (! $wb) or (! $yang_dir) ) {
    die "Format: ./wb_vs_yang_contigs.pl",
        " --wb|-w [WormBase assembly]",
        " --yang|-y [DIRECTORY of Yang contig files]",
        " --prefix|-p [prefix for Yang names]",
        "\n",
        ;
}

# Read assembly contigs; count nt in each non-N block; list in Yang-style:

open $WB, '<', $wb 
    or die "Can't open WormBase contigs file $wb: $!";

my $contig       = q{};
my $subcontig    = q{};
my $seq_index    = 0;
my $waiting      = 1;
my $nt_sizes_ref;

while (my $input = <$WB>) { 
    chomp $input;
    if ( $input =~ / \A > (\S+) /xms ) { 
        $contig = $1;
        $subcontig = q{};
        $seq_index = 0;
        $waiting   = 1;
    }
    if ( ( $input !~ / \A > /xms) 
         and ($input =~ / [acgtACGT]+ /xms ) ) { 
        $input =~ s/\s//g;
        my @residues = split //, $input;
        foreach my $res (@residues) { 
            # Test (! $waiting) condition *before* maybe changing $waiting value.
            if (! $waiting) {
                if ($res =~ /\A [acgtACGT] \z/xms ) {
                    $nt_sizes_ref->{'WB'}->{$subcontig}++;
                }
                elsif ($res =~ /\A [nN] \z/xms ) {
                    $waiting = 1;
                }
                else {
                    die "Can't parse residue $res in line $input\n";
                }
            }
            # And use elsif to further enforce either-or test.
            elsif ( ($waiting) and ($res =~ /\A [acgtACGT] \z/xms ) ) { 
                $waiting = 0;
                $seq_index++;
                $subcontig = $contig . '.' . $seq_index;
                $nt_sizes_ref->{'WB'}->{$subcontig}++;
            }
        }
    }
}

close $WB 
    or die "Can't close filehandle to WB contigs file $wb: $!";

# Read Shiaw-Pyng's multiple files from their directory; count nt in each subcontig; check vs. assembly.

opendir($YANG_DIR, $yang_dir) 
    or die "Can't opendir $yang_dir): $!";
while ( defined( $yang_file = readdir($YANG_DIR) ) ) { 
    if ( $yang_file !~ /\A \.{1,2} \z/xms ) { 
        my $yang_filepath = File::Spec->catfile($yang_dir,$yang_file);
        open $YANG_FILE, '<', $yang_filepath 
            or die "Can't open Shiaw-Pyng Yang contigs file $yang_file: $!";

        while (my $input = <$YANG_FILE>) { 
            chomp $input;
            if ( $input =~ / \A > (\S+) /xms ) { 
                $subcontig = $1;
                if ($prefix) { 
                    $subcontig = $prefix . $subcontig;
                }
            }
            elsif ( ( $input !~ / \A > /xms) 
                 and ($input =~ / \d+ /xms ) ) { 
                $input =~ s/\A\s+//;
                $input =~ s/\s+\z//;
                my @data_nos = split /\s+/, $input;
                my $res_count = @data_nos;
                $nt_sizes_ref->{'Yang'}->{$subcontig} += $res_count;
            }
        }
        close $YANG_FILE 
            or die "Can't close filehandle to Shiaw-Pyng",
                   " Yang contigs file $yang_file: $!",
                   ;
    }
}
closedir($YANG_DIR);

# Summarize all matches; report any discrepancies:

my $both_nt_equal    = 0;
my $wb_nt_unequal    = 0;
my $yang_nt_unequal  = 0;
my $wb_nt_only       = 0;
my $yang_nt_only     = 0;

my @output_lines = ();

foreach my $wb_subcontig ( sort keys %{ $nt_sizes_ref->{'WB'} } ) { 
    if (! exists $nt_sizes_ref->{'Yang'}->{$wb_subcontig} ) { 
        $wb_nt_only     += $nt_sizes_ref->{'WB'}->{$wb_subcontig};
        my $wb_size     = commify( $nt_sizes_ref->{'WB'}->{$wb_subcontig} );
        my $output_line = "$wb_subcontig"
                          . "\t"
                          . "$wb_size"
                          . "\t\t"
                          . "No Yang match"
                          . "\n"
                          ;
        push @output_lines, $output_line;
    }
    if ( exists $nt_sizes_ref->{'Yang'}->{$wb_subcontig} ) {
        my $wb_size     = commify( $nt_sizes_ref->{'WB'}->{$wb_subcontig} );
        my $yang_size   = commify( $nt_sizes_ref->{'Yang'}->{$wb_subcontig} );
        my $output_line = "$wb_subcontig"
                          . "\t"
                          . "$wb_size"
                          . "\t"
                          . "$yang_size"
                          . "\t"
                          ;
        if ( $wb_size ne $yang_size ) { 
            $output_line     .= "Diff. sizes\n";
            $wb_nt_unequal   += $nt_sizes_ref->{'WB'}->{$wb_subcontig};
            $yang_nt_unequal += $nt_sizes_ref->{'Yang'}->{$wb_subcontig};
        } 
        if ( $wb_size eq $yang_size ) { 
            $output_line .= "Ident. sizes\n";
            $both_nt_equal += $nt_sizes_ref->{'WB'}->{$wb_subcontig};
        }
        push @output_lines, $output_line;
    }
}

foreach my $yang_subcontig ( sort keys %{ $nt_sizes_ref->{'Yang'} } ) {
    if (! exists $nt_sizes_ref->{'WB'}->{$yang_subcontig}) {
        $yang_nt_only += $nt_sizes_ref->{'Yang'}->{$yang_subcontig};
        my $yang_size = commify($nt_sizes_ref->{'Yang'}->{$yang_subcontig});
        my $output_line = "$yang_subcontig"
                          . "\t\t"
                          . "$yang_size"
                          . "\t"
                          . "No WB match"
                          . "\n"
                          ;
        push @output_lines, $output_line;
    }
}

@output_lines = sort @output_lines;

$both_nt_equal    = commify($both_nt_equal);
$wb_nt_unequal    = commify($wb_nt_unequal);
$yang_nt_unequal  = commify($yang_nt_unequal);
$wb_nt_only       = commify($wb_nt_only);
$yang_nt_only     = commify($yang_nt_only);

print "\n";
print "Equally-sized subcontigs:       $both_nt_equal nt\n";
print "WB contigs, unequally sized:    $wb_nt_unequal nt\n";
print "Yang contigs, unequally sized:  $yang_nt_unequal nt\n";
print "WB contigs, no Yang match:      $wb_nt_only nt\n";
print "Yang contigs, no WB match:      $yang_nt_only nt\n";
print "\n";
print "Subcontig:\tWB (nt):\tYang (nt):\tStatus\n";

# Print unquoted to avoid automatic spacing between lines:
print @output_lines;
print "\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

