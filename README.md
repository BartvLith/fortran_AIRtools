# fortran_AIRtools
Partial translation of AIRtools II to fortran.

Algebraic Iterative Reconstruction (AIR) methods are a type of iterative linear solver that are popular in computed tomography (CT) reconstructions. Here's a breakdown of the various modules and the main program.

main.f90 is a program that implements a AIR methods. Using various flags and input files allows many different versions of Kaczmarz, random Kaczmarz, symmetric Kaczmarz and gradient descent to be used. Two stopping rules are also implemented: error gauge and mutual step. Use ./main.exe -h for help.

sparse_matrices.f90 implements Compressed Sparse Row (CSR) and Compressed Sparse Column (CSC) matrices together with various basic operations such as matrix-vector multiplication, transposition, etc.

system_generation.f90 currently implements some aspects of generating the CT linear system, such as generating the system matrix and reading in data. At the moment, this is limited to a fan linear geometry, but I intend to expand this in the future.

phantomgallery.f90 is intended to be a library of test problems for CT applications. Currently, there's only the Shepp-Logan phantom implemented, but I will add some more soon.

AIRMs.f90 contains the reconstruction methods. Currently, there is the Algebraic Reconstruction Technique family which contains Kaczmarz, symmetric Kaczmarz and random Kaczmarz. There's also an implementation of gradient descent. More methods may be added in the future, if I have time.