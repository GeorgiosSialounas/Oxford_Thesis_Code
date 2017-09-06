# Oxford_Thesis_Code

This repository contains the files from my thesis.  Not all this code is mine:

Code from Prof. P. Farrell:
- deflation.py
- petscsnessolver.py
- mynewton.py
- utils.py
- middelkamp2010.py

Code from Dr L. Liu:
- MSPIN_codes.m
- ls.m

My code:
everything else

The MSPIN code in python is middel_gs_mspinalg.py. 
You can ran both newton and MPSIN on this just comment out the one you don't want in lines 594-595.
The variational form and the Functions are given by the function in forward_mspin_total2.  The MSPIN
solver itself is mspin_alg_sep.  The newton solver finds 6 solutions with deflation. The MSPIN solver
finds only one.

For the matlabe examples I use mspin without linesearch on 5 algebraic problems in R^2.  The current file
numbering does not follow the ordering in my thesis.  I will correct these issue and upload a consistent version
as soon as I get access to matlab.  I will also add comments to the files.
