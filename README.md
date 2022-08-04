# Study of numerical methods for ODEs on the Hénon-Heiles differential equation system

I have implemented several numerical methods for solving ODEs, and have applied them to the Hénon-Heiles differential equation system to observe their properties. The methods are:

* Kutta's method of order 4
* Shampine-Bogacki method of order 3
* Kahan's method of order 2
* Störmer-Verlet method of order 2

The Hénon-Heiles system consists of the following set of equations:

$$
\begin{align*}
q_1' &= p_1 \\
q_2' &= p_2 \\
p_1' &= -q_1(1+2q_2) \\
p_2' &= -(q_2 +q_1^2 -q_2^2)
\end{align*}
$$

## Energy drift in non-symplectic numerical ODE solvers

By applying the following Hamilton function:

$$
\begin{align*}
H(q,p) &= \frac{1}{2}(p_1^2 +p_2^2) + U(q) \\
U(q) &= \frac{1}{2}(q_1^2 +q_2^2) + \lambda(q_1^2q_2-\frac{1}{3}q_2^3)
\end{align*}
$$

($\lambda = 1$, for ease) to our system of differential equations we can determine how well each of the implemented numerical methods can preserve the energy in the system. We will notice that only the symplectic methods (Kahans method and Störmer-Verlet) has this property, while our "traditonal" (non-symplectic) Runge-Kutta methods lacks the useful property.

## Poincaré mapping of the Hénon-Heiles system
Poncaré-mapping is often used to analyze dynamic, higher dimension, systems in simpler terms. The first main reason for calculating the Poincaré-map of the Hénon-Heiles system, was to determine if there existed a third invariant for the stellar motion inside a specific gravitational potential of a galaxy.

In our case we are going to map all the times $q_1$ in our system equals 0, at the same time as $p_1 > 0$. The method for calculating the Poncaré-map is to iterate through our system and notice when $q_1 = 0$ and $p_1 > 0$. Then we will use interpolation to calculate more accurate values of $q_2$ and $p_2$ which we then will scatter plot in a 2D-coordinate system, thus making it easier to analyze our 4-dimensional system.


The full report (with relatively messy python code(^:) can be found [here](https://github.com/arurans/MTFYMA/blob/master/TMA4320/Project%203/Project3%20without%20output.ipynb).

# Dependencies & Structure
This project was mostly done for me to try out numerical computation in C++ with the Armadillo and Eigen libraries, in addition to parallelization with OpenMP. All the computation is done in C++, while only the plotting part is done in Julia. This project uses the numerical libraries Armadillo and Eigen (separately) which must be installed, and OpenMP for parallelization (which should be a feature of new versions of gcc/g++ compilers, I used 9.4.0 for this project).


A setup tutorial for future me, if that time ever comes (Ubuntu 20.04):

To install armadillo (only openblas and lapack is required, but the rest is recommended):
```
sudo apt-get update && sudo apt-get upgrade
sudo apt install cmake libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev
sudo apt install libarmadillo-dev
```

To install Eigen:
```
sudo apt install libeigen3-dev
```

To set up the build folder (write armadillo instead of eigen if you prefer that):
```
git clone https://github.com/arurans/Henon-Heiles-system.git
cd eigen
mkdir build && cd "$_"
cmake -S ../ -B .
```
And the following to compile and run (from the build folder, either in armadillo or eigen):

```
make
./hhp
```

The document structure is explained below (same for both armadillo and eigen):

```
plots---------------------------------------------------- All the plots
src
|-- methods---------------------------------------------- Implemented numerical methods
|   |-- CMakeLists.txt
|   |-- kahans.cpp
|   |-- kahans.h
|   |-- rk4.cpp
|   |-- rk4.h
|   |-- sb.cpp
|   |-- sb.h
|   |-- sv.cpp
|   `-- sv.h
|-- problems--------------------------------------------- Computing functions
|   |-- CMakeLists.txt
|   |-- compute.cpp
|   |-- compute.h
|   |-- hamiltonian.cpp
|   |-- hamiltonian.h
|   |-- poincare.cpp
|   `-- poincare.h
|-- CMakeLists.txt
|-- constants.h------------------------------------------ Constants used for computation
|-- storage_info.cpp------------------------------------- File names to store computed data
|-- storage_info.h--------------------------------------- and the SKIP_STORAGE variable
|-- utils.cpp
`-- utils.h
CMakeLists.txt
main.cpp------------------------------------------------ main file/call desired functions
```

# Practical information

Eigen is both faster to compile (as it lets you choose which submodules to include), and also runs faster than armadillo. I have optimized Eigen a bit, and given it the "SKIP_STORAGE" property, which is just a variable that lets you store only every "SKIP_STORAGE-th" iteration to the final matrix, thus allowing for computation over larger intervals, even with limited RAM. SKIP_STORAGE shuold probably be set to 1 when computing the Poincaré maps. This variable can be found in the following file:

```
./eigen/src/storage_info.cpp
```

In addition I have noticed that for this particular code, both the armadillo and eigen implementations runs quite a bit faster when compiled with -O1 optimization flag rather than -O2 or -O3. In

```
./eigen/CMakeLists.txt
```

```
set(CMAKE_CXX_FLAGS -O1)
```