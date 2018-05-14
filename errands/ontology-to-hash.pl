#!/usr/bin/perl -w

# ontology_to_hash.pl
# 
# Erich Schwarz, 1/5/01
#
# Purpose: convert *.ontology tables from geneontology.org into
#   hash tables usable in automatic gene annotation.
# .............................................................
# 
# Input: like this, from component.ontology:
#
#    $Gene_Ontology ; GO:0003673
#     $cellular_component ; GO:0005575
#      %membrane ; GO:0016020
#       <integral membrane protein ; GO:0016021
#      %cell fraction ; GO:0000267
#       %insoluble fraction ; GO:0005626
#       %membrane fraction ; GO:0005624
#       %soluble fraction ; GO:0005625
# .............................................................
# 
# Output: like this --
# 
# %component function (and %function, and %process --
#      these names will all come from X.ontology file name)
# 
# $component(0003673) = Gene_Ontology
# ...
# $component(0016021) = integral membrane protein
# ...
# $component(0005625) = soluble fraction
# .............................................................

print "\nThis script is just beginning construction.\n\n";
