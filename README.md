# TODO: simulation set

The goal is to run an ensemble of simulations to "describe" the Hilbert space of action of supersonic ideal fluid (Euler equations).

* Let's run simulations with a different total energy, which can be achieved by varying left pressure $`pL \in (0.5, 1.0) (Sod shock tube uses `pL = 1.0` and `pR = 0.1`). Initial velocity and density remain untouched.

Second parameter to vary is the resolution of the simulation, which will lead to a scan over different viscosity strength.

## Suggestions

Let's shoot for
* 20 steps of $`pL \in (0.5, 1.0)`$
* 20 steps of number of cells $`\in (500, 2500)`$

which will lead to 400 datapoints.

Note, that our Hilbert space should be independent of the resolution of the simulation, so we will want to map (interpolate) the action fields on a common grid (2500 cells makes sense).

# Example

Execute the 1D ideal fluid shock tube simulation (generates 100 snapthots by default)

`python Euler1D_numpy.py 1000`

where the input parameter `1000` sets number of mesh cells.

Visualization of the output:

`python plot_Euler1D.py temp/output_00.txt`

`python plot_Euler1D.py temp/output_20.txt`

`python plot_Euler1D.py temp/output_40.txt`

`python plot_Euler1D.py temp/output_60.txt`

`python plot_Euler1D.py temp/output_80.txt`

`python plot_Euler1D.py temp/output_100.txt`

## Action update

$`\begin{equation}
\frac{d \mathcal{S}}{d t}(x, t) = \left(\partial_a x\right)(x, t) \rho(x, t) \left( \frac{1}{2}v(x, t)^2 - \varepsilon(x, t) \right),~(1)
\end{equation}`$

where specific internal energy

$`\begin{equation}
\varepsilon = \frac{p}{(\gamma - 1) \rho},~(2)
\end{equation}`$

and fluid coordinate evolves as

$`\begin{equation}
\partial_t a(x, t) = - v(x, t),~(3)
\end{equation}`$

and we approximate $`\partial_a x`$ in the lab frame as

$`\begin{equation}
\left( \partial_a x \right)(x_i) \approx \frac{x_{i+1} - x_i}{a(x_{i+1}) - a(x_i)}.~(4)
\end{equation}`$

## Ideal fluid example

Let's $`a`$ is the fluid coordinate, then the ideal fluid Lagrangian density function reads

$`\begin{equation}
\mathcal{L}^{if}(\rho_0(a), \dot{f}(a, \tau), \partial_a f(a, \tau)) = \rho_0(a) \left( \frac{1}{2}\dot{f}(a, \tau)^2 - \varepsilon\left( s_0(a), \frac{\rho_0(a)}{\partial_a f(a, \tau)} \right) \right),~(L5)
\end{equation}`$

where $`\varepsilon(s, \rho) = c \rho^{\gamma - 1} \exp(\alpha s)`$, where $`c`$, $`\gamma`$, and $`\alpha`$ are constants, is the internal energy potential depending on density $`\rho(a, t) = \frac{\rho_0(a)}{\partial_a f(a, \tau)}`$. Note that $`\rho_0(a)`$ (density in with respect to fluid coordinates) and $`s_0(a)`$ (isentropic processs) do not change with $`\tau`$.

The conjugate momentum of field $`f`$ defined by $`(L2)`$ from $`(L5)`$ is

$`\begin{equation}
\pi(a, \tau) = \partial_{\dot{f}} \mathcal{L}_f^{if}(a, \tau) = \rho_0(a) \dot{f}(a, \tau).~(L6) 
\end{equation}`$

Finally, we obtain the ideal fluid action density $`\mathcal{S}^{if}_f(a, \tau)`$ for a given function $`f`$ from $`(L3)`$ using ideal fluid Lagrangian density $`(L5)`$ by solving an ordinary differential equation

$`\begin{equation}
\frac{d \mathcal{S}_f^{if}}{d \tau}(a, \tau) = \rho_0(a) \left( \frac{1}{2}\dot{f}(a, \tau)^2 - \varepsilon\left( \frac{\rho_0(a)}{\partial_a f(a, \tau)} \right) \right).~(L7)
\end{equation}`$

**Remark:**
*It should be noted, that $`f(a, \tau)`$ is the evolving spatial position of a fluid element $`a`$ in time, $`f(a, \tau) = x(a, \tau)`$ in 1D*.

It can be shown that $`\dot{f}(a, \tau)`$ and $`\partial_a f(a, \tau)`$ are equivalent to Euler variables of ideal fluid

$`\begin{align}
\dot{f}(a, \tau) &= v(x(a, \tau), \tau),~(L8)
\\
\partial_a f(a, \tau) &= \frac{\rho_0(a)}{\rho(x(a, \tau), \tau)},~(L9)
\end{align}`$

where $`v(x, \tau)`$ and $`\rho(x, \tau)`$ are solution to Euler equations

$`\begin{align}
\rho \left( \partial_\tau v + v \partial_x v \right) &= - \partial_x p,
\\
\partial_\tau \rho + \partial_x \left( \rho v \right) &= 0,~(L10)
\\
\partial_\tau s + v \partial_x s &= 0,
\end{align}`$

with isentropic closure $`p = (\gamma - 1) \rho \varepsilon(s, \rho) = c_0 \rho^\gamma`$, $`c_0`$ constant in fluid coordinates.

Note that the realtion between Eulerian coordinate $`x`$ and Lagrangian (fluid) coordinate $`a`$ reads

$`\begin{align}
x(a, \tau) &= a + \int_0^\tau v(x(a, \tilde{\tau}), \tilde{\tau}) d\tilde{\tau},~(L11)
\\
a(x, \tau) &= x - \int_0^\tau v(x, \tilde{\tau}) d\tilde{\tau},~(L12)
\end{align}`$

where we assume initial condition $`x(a, 0) = a`$.

We conclude this ideal fluid example by rewritting action density from $`(L7)`$ into laboratory Eulerian coordinates by using $`(L12)`$ as

$`\begin{equation}
\frac{d \mathcal{S}_f^{if}}{d \tau}(x, \tau) = \rho_0(a(x, \tau)) \left( \frac{1}{2}v(x, \tau)^2 - \varepsilon(x, \tau) \right),~(L13)
\end{equation}`$

and conjugate momentum $`(L6)`$ in lab coordinates as

$`\begin{equation}
\pi(x, \tau) = \rho_0(a(x, t)) v(x, \tau).~(L14) 
\end{equation}`$

The action $`S^{if}[f](\tau)`$ can be obtained from the solution of $`(L7)`$ via $`(L4)`$.

Note, that, equivalently, can be obtained from solution of $`(L13)`$ via

$`\begin{equation}
S[f](\tilde{\tau}) = \int_\Omega \mathcal{S}_f^{if}(x, \tilde{\tau}) \frac{1}{J}d^3x,~(L15)
\end{equation}`$

where the integration is carried out over laboratory coordinates (not the fluid coordinates, requiring Jacobian scaling $`J = \partial_a x(a(x, \tau), \tau)`$).

# Euler1D Benchmark
A comparison of various programming languages solving a 1D hydrodynamics problem.

The programs implement a simple finite-difference solver (the [Lax-Friedrichs method](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method)) for the 1D [Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)) (inviscid, compressible hydrodynamics) in various languages to compare execution speed and syntax.

The benchmark is to solve the standard [Sod Shock Tube](https://en.wikipedia.org/wiki/Sod_shock_tube), with NX = 5000 points and a CFL parameter of 0.9. Disk output is only used initially to verify implementation correctness, but is suppressed for the benchmarks. Ten runs are carried out with each implementation to get an idea of the execution time variability (which is generally found to be small).

Current implementations are:

- C/C++
- Fortran 90 
- Java
- Python, using native nested lists and `for` loops
- Python, using Numpy arrays and vectorized operations
- Julia, native
- Julia, leveraging the [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) library
- Rust

The tests were executed on a personal computer with an Intel Core i7-9700F processor running Gentoo Linux with kernel 6.1.31. Compiler/interpreter versions were 12.3.1 for gcc (C/C++ and Fortran 90), OpenJDK 17.0.6 for Java, CPython 3.8.17 and 3.11.4 for Python (with Numpy 1.24.4), 1.8.5 for Julia and 1.69.1 for Rust.

Execution time is measured internally by each program by comparing the platform's wall clock or CPU clock at the start and end of the main program, and printed to the terminal, in seconds, as the sole program output. For Julia+LoopVec the code was run manually within the Julia REPL after executing `using LoopVectorization` and `include("Euler1D_opt.jl")`, so as to compile/optimize as much as possible before the tests are run, but the program still measures its own run time internally; the results of doing this are in line with those obtained with the `@benchmark` macro from [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).

Optimization flags were used where available (e.g. `O3` for gcc and `opt-level=3` for Rust); see the [`run_benchmarks.py`](https://github.com/meithan/Euler1D_Benchmark/blob/main/run_benchmarks.py) script for details. The [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) library was used for the Julia+LoopVec benchmark to vectorize/optimize the main loops (by simply prepending the `@turbo` macro); thanks to Luis Arcos ([LAlbertoA](https://github.com/LAlbertoA)) for the tip. The two Python benchmarks were executed with both Python 3.8 and 3.11; it turns out that 3.11 is substantially faster!

It was my first time writing Julia and Rust code so those implementations might be a bit rough. If you know how to further optimize any of them to make the comparison more fair, please let me know!

## Results

The results of the benchmark are presented in the following plots. The bar heights are the average execution time averaged over 10 runs for each case, normalized to the execution time of the fasest benchmark implementation (so far, Fortran). In the first plot the native Python implementations are shown on a separate scale.

#### **Linear** Y scale
![Lin scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_lin.png)


#### **Logarithmic** Y scale
![Log scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_log.png)

