# FAIR tools
High-performance computing implementation of (parts of) the Matlab package AIR Tools II. For more information on AIR Tools II, see the [research article](https://link.springer.com/article/10.1007/s11075-017-0430-x); the [AIR Tools II homepage](http://people.compute.dtu.dk/pcha/AIRtoolsII/index.html); or the [GitHub repository](https://github.com/jakobsj/AIRToolsII).

### To-do list
- Always add more test phantoms.
- Implement blocking support.
- Make matrix-free implementations.

### Added features
- Parallel matrix multiplications, SIRT methods and gradient descent are a bit faster now.
- Storage of the matrix is only feasible when a binary data dump is used. This makes the file itself unreadable by humans or any other computer program that doesn't know the contents beforehand.

## Algebraic Iterative Reconstruction methods
Algebraic Iterative Reconstruction (AIR) methods are a type of iterative linear solver that are popular in computed tomography (CT) reconstructions. The class includes Kaczmarz's method and other Algebraic Reconstruction Techniques (ART), Cimmino's method and other Simultaneous Iterative Reconstruction Techniques (SIRT).
Here's a breakdown of the various modules and the main program.

## Contents

### main.f90
main.f90 is a program that implements a AIR methods. Using various flags and input files allows many different versions of Kaczmarz, random Kaczmarz, symmetric Kaczmarz and gradient descent to be used. Two stopping rules are also implemented: error gauge and mutual step, they only work with the ART methods. Perhaps other stopping rules will be added in the future, but at the moment, this has lowest priority. Use ./main.exe -h for help.

### sparse_matrices.f90
sparse_matrices.f90 implements Compressed Sparse Row (CSR) and Compressed Sparse Column (CSC) matrices together with various basic operations such as matrix-vector multiplication, transposition, etc. There is a specialised algorithm for adding multiples of a single row to another vector with optional upper or lower bounds. This last one is a key component of Kaczmarz-type algorithms.

### system_generation.f90
system_generation.f90 currently implements some aspects of generating the CT linear system, such as generating the system matrix and reading in data. At the moment, this is limited to a fan linear geometry, but I intend to expand this in the future.

## AIRMs.f90
AIRMs.f90 contains the reconstruction methods. Currently, there is the Algebraic Reconstruction Technique (ART) family which contains Kaczmarz, symmetric Kaczmarz and random Kaczmarz. Some members of the Simultaneous Iterative Reconstruction Technique (SIRT) are implemented: Cimmino's method and SART (simultaneous ART). There's also an implementation of gradient descent. More methods may be added in the future.

To use a particular method, use the flag -m[methodname]. Standard is Kaczmarz. There are various options available, such as box constraints and relaxation parameter. To activate a lower bound, use -L[lbound]. To activate the upper bound, use -U[ubound]. The relaxation parameter is set using -R[relaxpar]. There are also options to use a stopping rule by means of -s[stoprulename] and ordering with -O[orderingfile]. I will provide some example programs  soon.

### Kaczmarz
Processes each row separately by projecting to the hyperplanes successively, methodname = 'kaczmarz'. Cyclicly goes through the row ordering. No input flag -m will also result in Kaczmarz being used with "natural" row ordering.

### Symmetric Kaczmarz
Each sweep of the rows is followed by a sweep in the reversed order. This has the effect of producing a symmetric iteration matrix. methodname = 'symkaczmarz' or 'sym'.

### Random Kaczmarz
The row ordering is picked randomly instead of cyclicly. There is an internal option to use row probability scaled with row norms, according to Strohmer and Vershynin. Otherwise uniform probability is used. methodname = 'randkaczmarz' or 'rand'.

### Cimmino
Processes the rows simultaneously by reflecting in the hyperplanes one by one and averaging the resulting intermediate solutions. methodname = 'cimmino'.

### SART
Processes the rows simultaneously by weighing the rows and columns to unit vectors. methodname = 'sart' or 'sirt'.

### Gradient descent
Classical gradient descend, computes the gradient of 1/2 \| A*x-b \|^2. methodname = 'graddescent' or 'gd'.


## phantomgallery.f90
phantomgallery.f90 is intended to be a library of test problems for CT applications. Currently there are some 9 phantoms implemented.
The flag -t[phantomname] will activate test mode in main and generate a phantom according to phantomname. It is possible to adjust the noise level with the flag -n[noiselevel].

### Shepp-Logan
The classic alien head phantom, phantomname = 'shepplogan'.

### Smooth
A smooth test phantom made by superposing a few Gaussians, phantomname = 'smooth'.

### Three phases
A phantom mimicking three phases of a liquid, phantomname = '3phases'.

### Smooth three phases
Three phase liquid with some diffusion, phantomname = '3phasesmooth'.

### Four phases
Four phase liquid, phantomname = '4phases'.

### Binary
Phantom mimicking stratified layers in a material, phantomname = 'binary'.

### Minimum Spanning Tree
A minimum spanning tree of a collection of Gaussian random points. The edges are drawn as squares and given a random value. phantomname = 'mst'.

### Grains
A collection of Voronoi cells mimicking with different values, phantomname = 'grains'.

### Bubbles
A collection of bubbles with different values, phantomname = 'bubbles'.



