# PolyMorph
A single-chain polymer conformation generator

PolyMorph is a C program with minimal dependencies (just the C-standard library) designed to build a single N-residue linear homopolymer chain defined from a single PDB (Protein Databank) formatted file and morph its coordinates to match arbitrary parametric curves.

NOTE: GPT 5.2 was used to assist in general code clean-up, removal of initial dependencies, and the implementation of the minimal RPN-parametric curve interpreter. 

How it works:
---------------------------------
1) Polymorph parses a given PDB and two selected atoms (to define the length of the monomer and identify its connection points)
2) A series of parametric functions are then defined and parsed using a simple RPN interpreter
3) The arc-length of the monomer is used to scale the defined parametric curve
4) The monomers are then placed along the curve.
