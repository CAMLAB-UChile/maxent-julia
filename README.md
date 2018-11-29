# maxent-julia
A Julia implementation of the maximum-entropy basis functions

IMPORTANT: Julia has recently updated to Julia 1.0. I noticed that maxent-julia is no longer compatible with Julia 1.0. Instead, use Julia 0.6.4 or Julia 0.7.0 to run maxent-julia. Note that using Julia 0.7.0 will throw lots of warnings. This is normal as these warnings are deprecation warnings that were intentionally implemented in Julia 0.7.0 to help Julia programmers to migrate to Julia v1.0. These warnings can be switched off as follows:

$julia --depwarn=no main.jl

# Author
<a href="https://github.com/aaortizb">Alejandro Ortiz-Bernardin</a>, Assistant Professor, Department of Mechanical Engineering, Universidad de Chile.

José M. Cáceres, Undergraduate Research Assistant, Department of Mechanical Engineering, Universidad de Chile.

# Instructions
The program is controlled by the main.jl function. This is the only function that
must be setup by the user. To execute the code, setup the problem parameters in
main.jl (further instructions are given there) and run it.

When setting up main.jl make sure that
  - length(x) = dim
  - if dim>=2, length(x) = size(ncoord)[2]
  - if dim=1, length(x) = length(ncoord)/n
  - size(ncoord)[1] = n
  - length(gamma) = n
  - length(ilambda) = dim

Anyway, an error is thrown when any of the previous equalities are not satisfied.

# License
This project is licensed under the GPL3 License. This program is free software; it can be redistributed or modified under the terms of the GNU General Public License 3 as published by the Free Software Foundation. 
