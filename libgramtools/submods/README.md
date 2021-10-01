## Submodules

This directory contains .cpp files that can each be compiled
into an executable.

They provide utility functionalities to gramtools.

* combine_jvcfs: merge jvcf JSONs into one
* encode_prg: convert a character-based description of a prg (e.g. A[T,C]G) into a 
  binary prg that can be used directly by gramtools (build)
* print_fm_index: from a character-based description of a prg, 
prints a table with its suffix array & BWT 
* visualise_prg: produce a graphviz dot graph file of a portion of a prg
