#!/usr/bin/env perl

# elements2fimohits.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/8/2011.
# Purpose: for each element, annotate FIMO hits of motifs (if any).

use strict;  
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use Roman;  # to get 'arabic' function

my %opts    = ();

# Make hash of motifs deemed nonredundant in Table 2 of Mortazavi et al. (2010):
my %nonredundant_set = ();
my @nonred_motifs           = qw( 1-1   1-2   1-3   1-4    1-5
                                  1-6         1-8   1-9    1-10
                                  1-11              1-14   1-15 
                                        1-17  1-18  1-19   1-20 
                                              2-18 
                                        2-22        2-24 
                                  2-26  2-27               2-30 );

foreach my $nonred_motif ( @nonred_motifs ) { 
    $nonredundant_set{$nonred_motif} = 1;
}

my $motif          = q{};
my $motif_set      = 0;
my $element        = q{};
my $q_value        = q{};
my $index_entity   = q{};
my $counted_entity = q{};

# Default is 1; must be set to be lower.
my $max_q_value = 1;

my $data_ref;

my $tab_elements;
my $tab_motifs;

my $compressed_tab;
my $only_best_q;

GetOptions ( 'compressed_tab'        => \$compressed_tab,
             'tab_elements|tab_eles' => \$tab_elements,
             'tab_motifs|tab_mots'   => \$tab_motifs,
             'max_q_value|max_q=f'   => \$max_q_value,
             'only_best_q'           => \$only_best_q,
             'motif1|m1|1=s'         => \$opts{'motif1'},
             'motif2|m2|2=s'         => \$opts{'motif2'}, 
             'prefix=s'              => \$opts{'prefix'},
             'help'                  => \$opts{'help'},     );

if (    $opts{'help'} 
     or (    ( (! $tab_elements ) and (! $tab_motifs ) ) 
          or ( $tab_elements      and $tab_motifs      ) ) 
     or (! $opts{'motif1'}                               ) 
     or (! $opts{'motif2'}                               ) ) { 
    die "\n",
        "Format: annotate_elements.pl\n",
        "        --compressed_tab|-c        [optional flattening of tabular output]\n",
        "        --tab_elements|--tab_eles\n",
        "   [or] --tab_motifs|--tab_mots    [either element- or motif-centric table]\n",
        "        --max_q_value|--max_q      [maximum allowed q-value]\n",
        "        --only_best_q|-o           [only list best q-value for a list of eles./mots. associated with mots./eles.; works best w/o --max_q arg.]\n",
        "        --motif1|--m1|-1           [FIMO results for motif series 1]\n",
        "        --motif2|--m2|-2           [FIMO results for motif series 2]\n",
        "        --prefix|-p                [prefix for elements in the FIMO reports; defaults to 'ce6_ct_v6step13_']\n",
        "        --help|-h                  [evoke this help message]\n",
        "\n",
        ;
}

# Default value for $opts{'prefix'} if no value supplied:
$opts{'prefix'} ||= 'ce6_ct_v6step13_';

# Reject non-numerical or absurdly valued maximum q-values:

if ( (! looks_like_number($max_q_value) ) or ( $max_q_value < 0 ) or ( $max_q_value > 1 ) ) { 
    die "Maximum q-value cannot be $max_q_value; must be defined floating-point number between 0 and 1!\n";
}

# given FIMO outputs from PS1010 paper, give each element its list (if any) of motif hits, ordered by q-value (smallest best and first).

# Sort hits by ascending q-value. 
# Print out name, tab, plaintext hit list, plaintext hit list with parenthetical statistics.

# Read each file.  Slurp data for wanted motifs; this needs a specific list for each motif file!

foreach my $motif_set_file ( $opts{'motif1'}, $opts{'motif2'} ) { 
    $motif_set++;
    open my $MOTIF_DATA, '<', $motif_set_file or die "Can't open FIMO results for motif series $motif_set, $motif_set_file!\n";
    while (my $input = <$MOTIF_DATA>) { 
        chomp $input;
        if ( $input !~ /\A \# /xms ) { 
            # Map the motifs to element names.
            if ( $input =~ /\A (\d+) \t $opts{'prefix'} (\S+) \t .+ \t (\S+) \t [ACGT]+ \z /xms ) { 
                $motif   = $1;
                $element = $2;
                $q_value = $3;
                $motif   = $motif_set . '-' . $motif;
                if ( $nonredundant_set{$motif} ) { 
                    $data_ref->{'element'}->{$element}->{$q_value}->{$motif} = 1;
                    $data_ref->{'motif'}->{$motif}->{$q_value}->{$element}   = 1;
                }
            }
        }
    }
    close $MOTIF_DATA or die "Can't close filehandle to FIMO results for motif series $motif_set, $motif_set_file!\n";
}


if ($tab_elements) {
    $index_entity   = 'element';
    $counted_entity = 'motif';
}
if ($tab_motifs) { 
    $index_entity   = 'motif';
    $counted_entity = 'element';
}

foreach my $entity1 ( sort { namevalue($a) <=> namevalue($b) } keys %{ $data_ref->{$index_entity} } ) { 
    my @q_values1 = grep { $_ <= $max_q_value } sort { $a <=> $b } keys %{ $data_ref->{$index_entity}->{$entity1} };
    if ($only_best_q) {
        @q_values1 = ( $q_values1[0] );
    }
    if (! $compressed_tab ) {
        foreach my $q_value1 (@q_values1) {
            foreach my $entity2 ( sort { namevalue($a) <=> namevalue($b) } keys %{ $data_ref->{'element'}->{$entity1}->{$q_value1} } ) {
                print "$entity1\t$q_value1\t$entity2\n";
            }
        }
    }
    if ($compressed_tab) { 
        my @annots = ();
        my $total_count = 0;
        my $count_word = q{};
        foreach my $q_value1 (@q_values1) { 
            my @entities2 = sort { namevalue($a) <=> namevalue($b) } keys %{ $data_ref->{$index_entity}->{$entity1}->{$q_value1} };
            my $entity2_count =  @entities2;
            $total_count += $entity2_count;
            $count_word = opt_pluralize($counted_entity, $entity2_count);
            $entity2_count = "$entity2_count $count_word";
            my $entity2_list = join ', ', @entities2;
            $entity2_list = "[q = $q_value1, $entity2_count: $entity2_list]";
            push @annots, $entity2_list;
        }
        $count_word = opt_pluralize($counted_entity, $total_count);
        my $annot_text = join '; ', @annots;
        $annot_text = "$total_count $count_word: " . $annot_text;
        # Don't bother printing out empty entries:
        if ($total_count) { 
            print "$entity1\t$annot_text\n";
        }
    }
}

sub namevalue { 
    # Note that this presupposes that no $_name1_2 will ever be bigger than 999_999; if *that* assumption fails [!], so does the script.
    my $_mot1   = $_[0];
    my $_mot1_1 = q{};
    my $_mot1_2 = q{};
    my $_sortvalue = 0;
    if ( $_mot1 !~ /\A (?: \d+\- | chr(?:I|II|III|IV|V|X)\. ) \d+ \z /xms ) {
        die "Can't parse motif name $_mot1 in motif-sorting!\n";
    }
    if ($_mot1 =~ /\A ( \d+\- | chr(?:I|II|III|IV|V|X)\. )  (\d+) \z /xms ) {
        $_mot1_1 = $1;
        $_mot1_2 = $2;
    }
    if ( $_mot1_1 =~ /\A chr(I|II|III|IV|V|X) \. \z /xms ) { 
        $_mot1_1 = $1;
        $_mot1_1 = arabic($_mot1_1);
    }
    if ( $_mot1_1 =~ /\A (\d+) \- \z /xms ) { 
        $_mot1_1 = $1;
    }
    if ( $_mot1_2 > 999_999 ) { 
        die "Can't compare motif name $_mot1 in motif-sorting, because its second term ($_mot1_2) is greater than 999,999!\n";
    }
    $_mot1_1 = $_mot1_1 * 1_000_000;
    $_sortvalue = $_mot1_1 + $_mot1_2;
    return $_sortvalue;
}

sub opt_pluralize {
    my ($_counted_entity, $_entity_count) = @_;
    my $_count_word = $_counted_entity ;
    if ( $_entity_count > 1 ) {
        $_count_word = $_count_word . 's';
    }
    return $_count_word;
}

