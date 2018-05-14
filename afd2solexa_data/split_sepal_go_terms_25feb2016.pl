#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $input_data       = q{};
my $best_terms       = q{};
my $next_best_terms  = q{};
my $only_batch_terms = q{};
my @deg_sets         = ();

my $header = "Gene\tGO_term_set_A\tGO_term_set_B\tGO_terms_batch_only\tGO_terms_other";

my $data_ref;

my $help;

GetOptions ( 'input_data=s'       => \$input_data,
             'best_terms=s'       => \$best_terms,
             'next_best_terms=s'  => \$next_best_terms,
             'only_batch_terms=s' => \$only_batch_terms,
             'deg_sets=s{,}'      => \@deg_sets,
             'help'               => \$help,        );

if ($help or (! $input_data) or (! $best_terms) or (! $next_best_terms) or (! $only_batch_terms) or (! @deg_sets) ) { 
    die "Format: split_sepal_go_terms_24feb2016.pl\n",
        "            -i|--input_data        [simple gene-to-GO term table]\n",
        "            -b|--best_terms        [GO terms w/ highest sig. for genotype changes]\n",
        "            -n|--next_best_terms   [GO terms w/ lower sig. for genotype changes]\n",
        "            -o|--only_batch_terms  [GO terms w/ only sig. for batch changes]\n",
        "            -d|--deg_sets          [one or more lists of genes sig. assoc. w/ changes]\n",
        "            -h|--help              [print this message]\n",
        "            [print complex annotation table to <STDOUT>]\n",
        ;
}

open my $BEST, '<', $best_terms;
while (my $input = <$BEST>) {
    chomp $input;
    parse_go_desc_file($input, 'best', $best_terms);
}
close $BEST;

open my $NEXT, '<', $next_best_terms;
while (my $input = <$NEXT>) {
    chomp $input;
    parse_go_desc_file($input, 'next', $next_best_terms);
}
close $NEXT;

open my $BATCH, '<', $only_batch_terms;
while (my $input = <$BATCH>) {   
    chomp $input;
    parse_go_desc_file($input, 'batch', $only_batch_terms);
}
close $BATCH;

# Get genes that we denoted as *significantly* changed in an up or down condition.
foreach my $deg_set (@deg_sets) {
    # Get a plain-text description of what the DEG set is about.
    my $deg_desc = basename($deg_set);

    # Revise the plain-text desc.
    if ( $deg_desc =~ /\A (\S+reg) \.gene_list\.txt \z/xms ) { 
        $deg_desc = $1;
        $deg_desc =~ s/ATML1__LGO/ATML1::LGO/g; 
        $deg_desc =~ s/_/ /g;

        # Mostly, the above is fine, but one term needs to *keep* '_':
        $deg_desc =~ s/Col WT/Col_WT/g;

        $deg_desc =~ s/\.vs\./ vs. /g; 
        $deg_desc =~ s/\.upreg/, upreg./g; 
        $deg_desc =~ s/\.downreg/, downreg./g; 
    }
    else {
        die "Cannot extract DEG type from name of $deg_set\n";
    }

    # Check mappability of the plain-text description to GO term significance conditions.
    # Do not make this lethal, because legitimate situations exist where this will fail! (e.g., very short gene lists).
    # But, do make errors easy to catch by enumerating those classes for which GO terms have been seen,
    #     so that trivial errors of different naming become obvious.

    if (! exists $data_ref->{'change_seen'}->{$deg_desc} ) {
        my @changes_seen = sort keys %{$data_ref->{'change_seen'} };
        warn "For $deg_set, cannot identify GO terms associated with putative change:\n";
        warn "    $deg_desc\n";
        warn "But have seen the following changes:\n";
        foreach my $change_seen (@changes_seen) {
            warn "    $change_seen\n";
        }
    }

    # Then, link genes in that set to the description.
    open my $DEG, '<', $deg_set;
    while (my $input = <$DEG>) {
        chomp $input;
        if ( $input =~ / \A (AT (?: 1|2|3|4|5|C|M) G\d+) \z/xms ) {
            my $tair_gene = $1;
            $data_ref->{'tair_gene'}->{$tair_gene}->{'deg_desc'}->{$deg_desc} = 1;
        }
        else {
            die "From DEG set list $deg_set, cannot parse gene name $input\n";
        }
    }
    close $DEG;
}

open my $DATA, '<', $input_data;
while (my $input = <$DATA>) {
    chomp $input;
    if ( $input =~ / \A (AT (?: 1|2|3|4|5|C|M) G\d+) \t ([^\t]+) \z/xms ) { 
        my $tair_gene = $1;
        my $go_text   = $2;
        my @go_terms  = split /; /, $go_text;

        my @best_go_terms  = ();
        my @next_go_terms  = ();
        my @batch_go_terms = ();
        my @other_go_terms = ();

        foreach my $go_term (@go_terms) {
            if ( exists $data_ref->{'best'}->{$go_term} ) {
                if ( exists $data_ref->{'next'}->{$go_text} ) {
                    die "GO term $go_term annotated both as 'best' and as 'next'\n";
                }
                if ( exists $data_ref->{'batch'}->{$go_term} ) {
                    die "GO term $go_term annotated both as 'best' and as 'batch'\n";
                }
                push @best_go_terms, $go_term;
            }

            elsif ( exists $data_ref->{'next'}->{$go_term} ) {
                if ( exists $data_ref->{'batch'}->{$go_term} ) {
                    die "GO term $go_term annotated both as 'next' and as 'batch'\n";
                }
                push @next_go_terms, $go_term;
            }

            elsif ( exists $data_ref->{'batch'}->{$go_text} ) {
                push @batch_go_terms, $go_term;
            }

            else { 
                push @other_go_terms, $go_term;
            }
        }            

        @best_go_terms  = @{ sort_uniq_noblank(\@best_go_terms) };
        @next_go_terms  = @{ sort_uniq_noblank(\@next_go_terms) };
        @batch_go_terms = @{ sort_uniq_noblank(\@batch_go_terms) };
        @other_go_terms = @{ sort_uniq_noblank(\@other_go_terms) };

        my @best_go_texts  = ();
        my @next_go_texts  = ();
        my @batch_go_texts = ();

        go_term_list2_go_text_list($tair_gene, \@best_go_terms, 'best', \@best_go_texts);
        go_term_list2_go_text_list($tair_gene, \@next_go_terms, 'next', \@next_go_texts);
        go_term_list2_go_text_list($tair_gene, \@batch_go_terms, 'batch', \@batch_go_texts);

        my $best_go_text  = join '; ', @best_go_texts;
        my $next_go_text  = join '; ', @next_go_texts;
        my $batch_go_text = join '; ', @batch_go_texts;
        my $other_go_text = join '; ', @other_go_terms;

        print "$header\n" if $header;
        $header = q{};

        print "$tair_gene\t$best_go_text\t$next_go_text\t$batch_go_text\t$other_go_text\n";
    }
    elsif ( $input !~ /\A Gene \t GO_term \z/xms ) {
        die "From input data file $input_data, cannot parse: $input\n";
    }
}
close $DATA;

sub parse_go_desc_file {
    my $_input  = $_[0];
    my $_type   = $_[1];
    my $_source = $_[2];

    chomp $_input;

    if ( $_input !~ /\A GO[ _]term \t p-value\[condition\] \z/xms ) {
        if ( $_input =~ / \A ([^\t]+) \t ([^\t]+) \z/xms ) {
            my $_go_term   = $1;
            my $_cond_text = $2;
            if ( exists $data_ref->{$_type}->{$_go_term}->{'change'} ) {
                die "In $_type GO term file $_source, redundant conditions associated with $_go_term\n";
            }
        
            my @_conds = split /; /, $_cond_text;
            foreach my $_cond (@_conds) {
                $_cond =~ s/\A\s+//;
                $_cond =~ s/\s+\z//;
            
                my $_change  = q{};
                my $_p_value = q{};
                if ( $_cond =~ /\A (\S+) \s+ \[ ([^\[\]]+) \] \z/xms ) {
                     $_p_value = $1;
                     $_change  = $2;
                }
                else {
                    die "In $_type GO term file $_source, cannot parse condition: $_cond\n";
                }
            
                $data_ref->{$_type}->{$_go_term}->{'change'}->{$_change}->{'_p_value'} = $_p_value;
                $data_ref->{'change_seen'}->{$_change} = 1;
            }
        }
        else {
            die "In $_type GO term file $_source, cannot parse: $_input\n";
        }
    }
    return;
}

sub sort_uniq_noblank {
    my @_orig_array = @{ $_[0] };
    @_orig_array = grep { /\S/ } @_orig_array;
    @_orig_array = sort @_orig_array; 
    @_orig_array = uniq(@_orig_array);
    return \@_orig_array;
}

sub go_term_list2_go_text_list {
    my $_tair_gene         = $_[0];
    my @_list_go_terms     = @{ $_[1] };
    my $_type              = $_[2];
    my $_list_go_texts_ref = $_[3];

    foreach my $_go_term (@_list_go_terms) {
        my @_changes       = sort keys %{ $data_ref->{'tair_gene'}->{$_tair_gene}->{'deg_desc'} };
        my @_sig_changes   = ();
        my %_seen_sig      = ();
        my @_insig_changes = ();
        my @_sig_descs     = ();
        foreach my $_change (@_changes) {
            if ( exists $data_ref->{$_type}->{$_go_term}->{'change'}->{$_change}->{'p_value'} ) {
                push @_sig_changes, $_change;
            }
        }

        foreach my $_sig_change (@_sig_changes) {
            $_seen_sig{$_sig_change} = 1;
        }
        foreach my $_change (@_changes) {
            if (! $_seen_sig{$_change} ) {
                push @_insig_changes, $_change;
            }
        }
        @_insig_changes = sort @_insig_changes;
        @_insig_changes = uniq(@_insig_changes);

        @_sig_descs = map { "$data_ref->{$_type}->{$_go_term}->{'change'}->{$_}->{'p_value'} [$_]" }
                      sort {     $data_ref->{$_type}->{$_go_term}->{'change'}->{$a}->{'p_value'}
                             <=> $data_ref->{$_type}->{$_go_term}->{'change'}->{$a}->{'p_value'} }
                      @_sig_changes;
        my @_all_descs = ();
        push @_all_descs, @_sig_descs;
        push @_all_descs, @_insig_changes;

        my $_output_go_text = join ', ', @_all_descs;
        if ( $_output_go_text !~ /\S/xms ) {
            $_output_go_text = 'not sig. gene';
        }
        $_output_go_text = "$_go_term ($_output_go_text)";
        push @{ $_list_go_texts_ref }, $_output_go_text;
    }

    @{ $_list_go_texts_ref } = sort @{ $_list_go_texts_ref };
    @{ $_list_go_texts_ref } = uniq @{ $_list_go_texts_ref };
    return;
}

