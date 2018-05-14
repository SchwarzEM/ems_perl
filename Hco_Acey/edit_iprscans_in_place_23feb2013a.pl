#!/usr/bin/env perl

use strict;
use warnings;

while (my $main_script = <>) { 
    chomp $main_script;
    my $backup_script = $main_script . '.backup_22feb2013b';
    my $inline_perl_text_patch = q{ perl -ne ' s/mem\=10gb/mem=7gb/; print; ' };
    system qq{    mv -i $main_script $backup_script ;\n};
    system qq{    cat $backup_script | grep -v "feature=intel10" | $inline_perl_text_patch > $main_script ;\n};
}

