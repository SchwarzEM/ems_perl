#!/usr/bin/env perl

# strict_vs_exp_omcls.pl -- Erich Schwarz <ems394@cornell.edu>, 3/15/2013; minor tweak on 1/31/2017 to allow (OrthoMCL-reformatted) OrthoFinder input; another tweak on 3/23/2020 to allow Drosophila gene names like "l(3)wacko".
# Purpose: given OrthoMCL text (or OrthoMCL-like OrthoFinder text), select for taxa with a variety of conditions. Originally a script to 'allow one or more taxa to expand'.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $input_omcl = q{};

my @only      = ();
my %only_taxa = ();

my @required      = ();
my %required_taxa = ();

my @permitted      = ();
my %permitted_taxa = ();

my @banned      = ();
my %banned_taxa = ();

my @variable   = ();
my %var_taxa   = ();

my @expandable = ();
my %exp_taxa   = ();

my @multiple   = ();
my %mult_taxa  = ();

my $all_unspec_var;

# This is used to list all taxa *specified* in the arguments, which implies which taxa can have 'all_unspec_var' applied to them.
my %spec_taxa  = ();

my $help;

GetOptions ( 'input_omcl:s'     => \$input_omcl,
             'only=s{1,}'       => \@only,
             'required=s{1,}'   => \@required,
             'permitted=s{1,}'  => \@permitted,
             'banned=s{1,}'     => \@banned,
             'variable=s{1,}'   => \@variable,
             'expandable=s{1,}' => \@expandable,
             'multiple=s{1,}'   => \@multiple,
             'all_unspec_var'   => \$all_unspec_var,
             'help'             => \$help,       );

if ( ($help) or (! $input_omcl ) ) { 
     die "Format: strict_vs_exp_omcls.pl\n",
         "            --input_omcl|-i      [input: either OrthoMCL text file, or '-' for input stream]\n",

         "            --only|-o            [require that orthology groups have all of these species, and only them;",
         " only allows 1 instance/group by default, but this can be overriden by -e or -m]\n",

         "            --required|-r        [require that orthology groups have all of these species, but others allowed;",
         " only allows 1 instance/group by default, but this can be overriden by -e or -m]\n",

         "            --permitted|-p       [require that each orthology group has at least one of these species",
         " (but need not be the same 1+ spp. between groups); currently does not give option to enforce maximum of 1/group]\n",

         "            --banned|-b          [ban any orthology group which has any one or more of these species]\n",
         "            --variable|-v        [allow 0-N instances/group for these taxa -- i.e., completely variable]\n",
         "            --expandable|-e      [allow and require 1+ instances/group for these taxa]\n",
         "            --multiple|-m        [allow and require 2+ instances/group for these taxa]\n",
         "            --all_unspec_var|-a  [make all taxa not otherwise specified -v (shortcut command, useful with -r or -p)]\n",
         ;
}

foreach my $only_taxon (@only) {
    $only_taxa{$only_taxon} = 1;
    $spec_taxa{$only_taxon} = 1;
}

foreach my $required_taxon (@required) {
    if ( exists $only_taxa{$required_taxon} ) {
        die "Cannot simultaneously have 'only' and 'required' taxon $required_taxon\n";
    }
    $required_taxa{$required_taxon} = 1;
    $spec_taxa{$required_taxon} = 1;
}

foreach my $permitted_taxon (@permitted) {
    if (    ( exists $only_taxa{$permitted_taxon} ) 
         or ( exists $required_taxa{$permitted_taxon} ) ) {
        die "Cannot simultaneously require and permit taxon $permitted_taxon\n";
    }
    $permitted_taxa{$permitted_taxon} = 1;
    $spec_taxa{$permitted_taxon} = 1;
}

foreach my $banned_taxon (@banned) {
    if (    ( exists $only_taxa{$banned_taxon}      ) 
         or ( exists $required_taxa{$banned_taxon}  ) 
         or ( exists $permitted_taxa{$banned_taxon} ) ) {
        die "Cannot simultaneously ban and require/permit taxon $banned_taxon\n";
    }
    $banned_taxa{$banned_taxon} = 1;
    $spec_taxa{$banned_taxon} = 1;
}

foreach my $var_taxon (@variable) {
    if (    ( exists $only_taxa{$var_taxon}      ) 
         or ( exists $required_taxa{$var_taxon}  )
       	 or ( exists $permitted_taxa{$var_taxon} ) 
         or ( exists $banned_taxa{$var_taxon}    ) ) {
        die "Cannot allow 'variable' taxon $var_taxon while also requiring, permitting, or banning it\n";
    }
    $var_taxa{$var_taxon}  = 1;
    $spec_taxa{$var_taxon} = 1;
}

foreach my $exp_taxon (@expandable) {
    if (    ( exists $banned_taxa{$exp_taxon} ) 
         or ( exists $var_taxa{$exp_taxon}    ) ) { 
        die "Cannot have expandable taxon $exp_taxon that is also banned or variable\n";
    }
    $exp_taxa{$exp_taxon} = 1;
    $spec_taxa{$exp_taxon} = 1;
}

foreach my $mult_taxon (@multiple) {
    if (    ( exists $banned_taxa{$mult_taxon} ) 
         or ( exists $var_taxa{$mult_taxon}    ) 
         or ( exists $exp_taxa{$mult_taxon}    ) ) { 
        die "Cannot have multiple taxon $mult_taxon that is also banned, variable, or expandable\n";
    }
    $mult_taxa{$mult_taxon} = 1;
    $spec_taxa{$mult_taxon} = 1;
}

my $INPUT_OMCL;
if ($input_omcl eq '-') {
    # Special case: get the stdin handle
    $INPUT_OMCL = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_OMCL, '<', $input_omcl or die "Cannot open OrthoMCL file $input_omcl: $!";
}

while (my $input = <$INPUT_OMCL>) { 
    # Each line has one orthology group.  Decisions about whether to print them are done line by line.
    chomp $input;

    # The default choice *is* to print, but this choice can be irreversibly negated for many reasons.
    my $printable = 1;

    # For each orthology group, these start at zero:
    my %counted_taxon       = ();
    my @omcl_seen_taxa      = ();
    my @omcl_expanded_taxa  = ();
    my @omcl_multiple_taxa  = ();
    my @omcl_permitted_taxa = ();

    # Obligatory parsing of the input line.
    if ($input =~ / \A
                    [A-Z]+\d+                     # no longer demands 'ORTHOMCL' as only permissible starting text 
                    \( \d+\sgenes,\d+\staxa \)    # just (punctuation)
                    : \s+ \b                      # Enforce word boundary just before start of $oseqs_line.

                    (\S.+\S)                      # $1 -> $oseqs_line

                    \s* 
                    \z 
                  /xms) {

        # Note that we refer to '$oseqs_line', not '$oprots_line', 
        #     because the OrthoMCL line may have either proteins or genes.

        my $oseqs_line = $1;
        my @orthoseqs = split /\s+/, $oseqs_line;

        foreach my $o_seq (@orthoseqs) { 
            # Obligatory parsing of the $o_seq set; failure to parse dies.
 
            # Formerly used:
            # if ( $o_seq =~ / \A [^\s\(\)]+ \( ( [^\s\(\)]+ ) \) \z /xms ) {
            # But Drosophila has weird gene/protein names with internal parentheses, so...

            if ( $o_seq =~ / \A \S+ \( ( [^\s\(\)]+ ) \) \z /xms ) { 
                my $species = $1;

                # Always count the species, instance by instance through the OrthoMCL group.
                $counted_taxon{$species} += 1;

                # Assign the taxon to 'var' (and then mark it as 'spec'!), 
                #     if not otherwise specified and $all_unspec_var is in effect:
                if ($all_unspec_var and (! exists $spec_taxa{$species} ) ) {
                    $var_taxa{$species}  = 1;
                    $spec_taxa{$species} = 1;
                }

                # Censor the group if a species in that group is banned; 
                #     or if the species is not 'only', 'required', or 'permitted'; 
                #     or if it is not named in 'var'/'exp'/'mult'.
                if (    ( exists $banned_taxa{$species}              ) 
                     or (     (! exists $only_taxa{$species}      )
                          and (! exists $required_taxa{$species}  )
                          and (! exists $permitted_taxa{$species} )
                          and (! exists $var_taxa{$species}       )
                          and (! exists $exp_taxa{$species}       )
                          and (! exists $mult_taxa{$species}      )  ) 
                   ) {
                    $printable = 0;
                }
            }

            # Enforce correct format of 'Parse the $o_seq set':
            else {
                die "Can't parse sequence $o_seq in input $input!\n";
            }
        }
    }

    # Enforce correct format and success of 'parse the input line':
    else { 
        die "From OrthoMCL file $input_omcl, can't parse input line: $input\n";
    }

    @omcl_seen_taxa      = sort keys %counted_taxon;
    @omcl_expanded_taxa  = sort keys %exp_taxa;
    @omcl_multiple_taxa  = sort keys %mult_taxa;
    @omcl_permitted_taxa = sort keys %permitted_taxa ;

    # Require that expanded taxa have >=1 representatives in the orthology group:
    foreach my $omcl_expanded_taxon (@omcl_expanded_taxa) { 
        if ( (! $counted_taxon{$omcl_expanded_taxon} ) or ( $counted_taxon{$omcl_expanded_taxon} < 1 ) ) {
            $printable = 0;
        }
    }

    # Require that multiple taxa have >=2 representatives in the orthology group:
    foreach my $omcl_multiple_taxon (@omcl_multiple_taxa) {
        if ( (! $counted_taxon{$omcl_multiple_taxon} ) or ( $counted_taxon{$omcl_multiple_taxon} < 2 ) ) {
            $printable = 0;
        }
    }

    # Require that all 'only' or 'required' taxa be represented in the orthology group.
    my @obligatory_taxa = @only;
    push @obligatory_taxa, @required;

    # Prevent q{} from getting classed as a 'taxon', 
    #     which seriously screws up the program's parsing logic...
    @obligatory_taxa = grep { /\S/ } @obligatory_taxa;

    foreach my $obligatory_taxon (@obligatory_taxa) { 
        if (! $counted_taxon{$obligatory_taxon} ) {
            $printable = 0;
        }
    }

    # Require that 'only' and 'permitted' taxa have no more than 1 representative 
    #     in each orthology group, unless otherwise permitted:
    my @allowed_taxa = ();

    push @allowed_taxa, @only;
    push @allowed_taxa, @required;
    push @allowed_taxa, @permitted;

    # Again, death to q{} as a 'taxon'.
    @allowed_taxa = grep { /\S/ } @allowed_taxa;

    foreach my $allowed_taxon (@allowed_taxa) {
        if ( exists $counted_taxon{$allowed_taxon} ) {
            if ( $counted_taxon{$allowed_taxon} != 1 ) {  
                # It requires *all* of these three failures to disqualify an OrthoMCL group; any one success leave it $printable = 1!
                if (     (! exists $exp_taxa{$allowed_taxon}       ) 
                     and (! exists $mult_taxa{$allowed_taxon}      ) 
                     and (! exists $permitted_taxa{$allowed_taxon} )
                     and (! exists $var_taxa{$allowed_taxon}       ) ) { 
                    $printable = 0;
                }
            }
        }
    }

    # Require that, if there *are* 'only' taxa, that each observed taxon belong to that set, *and* vice versa.
    # I.e., an orthology group with only one species from 'only' should be disqualified!
    if (@only) { 
        foreach my $omcl_seen_taxon (@omcl_seen_taxa) {
            if (! exists $only_taxa{$omcl_seen_taxon} ) { 
                $printable = 0;
            }
        }
        foreach my $only_taxon (@only) { 
            if (! exists $counted_taxon{$only_taxon} ) { 
                $printable = 0;
            }
        }
        foreach my $required_taxon (@required) {
            if (! exists $counted_taxon{$required_taxon} ) {
                $printable = 0;
            }
        }
    }

    # If an orthology group has passed all of the above tests, then ...
    # Require that permitted taxa *collectively* have >1 representatives in the orthology group.
    if ($printable and @omcl_permitted_taxa) {
        # Abolish printable status for this orthology group, but restore it *if* at least one permittable taxon exists in the group.
        $printable = 0;
        foreach my $omcl_permitted_taxon (@omcl_permitted_taxa) {
            if ( $counted_taxon{$omcl_permitted_taxon} ) {
                $printable = 1;
            }
	}
    }

    # Finally, print any lines that passed *all* of the above tests.
    if ($printable) { 
        print "$input\n";
    }
}

close $INPUT_OMCL or die "Cannot close filehandle to OrthoMCL file $input_omcl: $!";

