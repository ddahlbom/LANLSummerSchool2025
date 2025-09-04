# LANLSummerSchool2025

This repo contains a number of very simple [Sunny](https://github.com/SunnySuite/Sunny.jl) scripts illustrating points mentioned in the lecture, "A Brief introduction to Classical and Semiclassical Spin Dynamics," given at the Computational Condensed Matter Summer School hosted by Los Alamos National Laboratory in June, 2025.

## Executing example
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> LANLSummerSchool2025

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "LANLSummerSchool2025"
```
which auto-activate the project and enable local path handling from DrWatson.
