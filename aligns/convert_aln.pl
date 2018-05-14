#!/usr/bin/env perl

# convert_aln.pl -- from Bio::AlignIO perldoc, by Erich Schwarz <emsch@its.caltech.edu>, 11/19/2011.
# Purpose: convert various sequence alignments; generalization of converting FASTA to Stockholm alignments, for (e.g.) HMMER.

use strict;
use warnings;
use Bio::AlignIO;
use File::Basename;
use Getopt::Long;

my @aligns     = q{};
my $in_format  = q{};
my $out_format = q{};
my $help;

my %aln2suffix = ( 'bl2seq'    => 'bl2seq',  # Bl2seq Blast output
                   'clustalw'  => 'aln',     # clustalw (.aln) format
                   'emboss'    => 'emboss',  # EMBOSS water and needle format
                   'fasta'     => 'fa',      # FASTA format
                   'maf'       => 'maf',     # Multiple Alignment Format
                   'mase'      => 'mase',    # mase (seaview) format
                   'mega'      => 'mega',    # MEGA format
                   'meme'      => 'meme',    # MEME format
                   'msf'       => 'msf',     # msf (GCG) format
                   'nexus'     => 'nexus',   # Swofford et al NEXUS format
                   'pfam'      => 'pfam',    # Pfam sequence alignment format
                   'phylip'    => 'phylip',  # Felsenstein PHYLIP format
                   'prodom'    => 'prodom',  # prodom (protein domain) format
                   'psi'       => 'psi',     # PSI-BLAST format
                   'selex'     => 'sel',     # selex (hmmer) format
                   'stockholm' => 'sto',     # stockholm format
);

GetOptions ( 'aligns=s{,}' => \@aligns,
             'in_form=s'   => \$in_format,
             'out_form=s'  => \$out_format,
             'help'        => \$help, );

if ($help or (! $in_format) or (! $out_format) or (! @aligns) ) { 
    print "Format: convert_aln.pl",
          " --align|-a [input alignment(s)]",
          " --in_form|-i [input format]",
          " --out_form|-o [output format]",
          " --help|-h [this message]",
          "\n",
          ;
    print_formats();
    exit;
}

if ( (! exists $aln2suffix{$in_format} ) or (! exists $aln2suffix{$out_format} ) ) { 
    print "Input format $in_format not accepted\n"   if (! exists $aln2suffix{$in_format} );
    print "Output format $out_format not accepted\n" if (! exists $aln2suffix{$out_format} );
    print_formats();
    exit;
}

foreach my $input_align (@aligns) { 
    my $output_align = basename($input_align);
    $output_align =~ s/\.$aln2suffix{$in_format}\z//;
    $output_align .= ".$aln2suffix{$out_format}";

    if (-r $input_align) { 
        if (-e $output_align ) { 
            die "Won't overwrite $output_align!\n";
        }
        my $in  = Bio::AlignIO->new( -file   => "$input_align",
                                     -format => "$in_format",        );
        my $out = Bio::AlignIO->new( -file   => ">$output_align", 
                                     -format => "$out_format",    );
        while ( my $aln = $in->next_aln() ) { 
            $out->write_aln($aln);
        }
    }
}

sub print_formats { 
    print "Supported formats, via Bio::AlignIO:\n";
    print "\n";
    print "    bl2seq:     Bl2seq Blast output\n";
    print "    clustalw:   clustalw (.aln) format\n";
    print "    emboss:     EMBOSS water and needle format\n";
    print "    fasta:      FASTA format\n";
    print "    maf:        Multiple Alignment Format\n";
    print "    mase:       mase (seaview) format\n";
    print "    mega:       MEGA format\n";
    print "    meme:       MEME format\n";
    print "    msf:        msf (GCG) format\n";
    print "    nexus:      Swofford et al NEXUS format\n";
    print "    pfam:       Pfam sequence alignment format\n";
    print "    phylip:     Felsenstein PHYLIP format\n";
    print "    prodom:     prodom (protein domain) format\n";
    print "    psi:        PSI-BLAST format\n";
    print "    selex:      selex (hmmer) format\n";
    print "    stockholm:  stockholm format\n";
    print "\n";
}
