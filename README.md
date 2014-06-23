alignment-algos
===============

C++98 [updated to C++11] implementation of my graduate work in protein sequence alignment tools

Summary
=======

This library produces a suite of tools that construct optimal and sets of near-optimal alignments based on minimizing edit distances between sequences using a variation on the Smith-Waterman method. It also contains tools that create reports on a number of geometric parameters for protein structures, such as surface area burial, hydrogen bonding residues, residue-to-residue distances, and more.

Features
========

  1. Near-optimal sequence alignment
  2. Template library design

Near-optimal sequence alignment
-------------------------------

A unique feature of the library is its ability to enumerate alignments rather than just optimizing them. This allows the protein structure modeler to visualize where sequence alignment variation can occur and focus attention on those regions.

Template library design
-----------------------

The core of the library is a dynamic programming algorithm that can be reused regardless of the sequence types and scoring function being used to evaluate the alignment. The core alignment algorithm can be swapped out itself!

Notes
=====

This work was developed in the Honig lab @ Columbia University, and some of the tools depend on components that I'm not sure have been made available to the public domain.
