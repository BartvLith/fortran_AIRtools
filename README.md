# FAIRtools
Partial translation of AIRtools II to fortran.

Algebraic Iterative Reconstruction (AIR) methods are a type of iterative linear solver that are popular in computed tomography (CT) reconstructions. Here's a breakdown of the various modules and the main program.

main.f90 is a program that implements a AIR methods. Using various flags and input files allows many different versions of Kaczmarz, random Kaczmarz, symmetric Kaczmarz and gradient descent to be used. Two stopping rules are also implemented: error gauge and mutual step, they only work with the ART methods. Perhaps other stopping rules will be added in the future, but at the moment, this has lowest priority. Use ./main.exe -h for help.

sparse_matrices.f90 implements Compressed Sparse Row (CSR) and Compressed Sparse Column (CSC) matrices together with various basic operations such as matrix-vector multiplication, transposition, etc. There is a specialised algorithm for adding multiples of a single row to another vector with optional upper or lower bounds. This last one is a key component of Kaczmarz-type algorithms.

system_generation.f90 currently implements some aspects of generating the CT linear system, such as generating the system matrix and reading in data. At the moment, this is limited to a fan linear geometry, but I intend to expand this in the future.

phantomgallery.f90 is intended to be a library of test problems for CT applications. There are a few phantoms in there. I will be adding more.

AIRMs.f90 contains the reconstruction methods. Currently, there is the Algebraic Reconstruction Technique (ART) family which contains Kaczmarz, symmetric Kaczmarz and random Kaczmarz. Some members of the Simultaneous Iterative Reconstruction Technique (SIRT) are implemented: Cimmino's method and SART (simultaneous ART). There's also an implementation of gradient descent. More methods may be added in the future.