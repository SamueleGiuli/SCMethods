# SCMethods (SCM)

A Fortran library providing self-consistency methods in one simple container.

## External Libraries
The SCM library depends upon the following libraries.
- [STDLIB](https://github.com/fortran-lang/stdlib)
- LAPACK
- BLAS
	
## Methods (to be improved)

At the present state the library give access to the following methods:

| method         | real32|real64|
|--------|-------|-------|
| linear_mixing  |   ✅  | ✅    | 
| DIIS           |   ✅  | ✅    |
	
## Building and using

The project was built using the [Fortran Package Manager](https://github.com/fortran-lang/fpm).

To use `SCMethods` within your FPM project, add the following to your `fpm.toml` file:
```toml
[dependencies]
SCMethods = { git="https://github.com/SamueleGiuli/SCMethods.git" }
```

⚠️ Please note that **currently the library cannot be built with the default compiler flags** requested by `fpm`.    
For GNU Fortran (the Fortran compiler of the GCC suite) you can use the 
 ```
 --flag "-ffree-line-length-none -fPIC -w -fallow-argument-mismatch -O3 -funroll-loops"
 ```
 option after `fpm build`, `fpm run`, `fpm test` or you can setup your shell to remember the choosen flags as
 ```
 export FPM_FFLAGS="-ffree-line-length-none -fPIC -w -fallow-argument-mismatch -O3 -funroll-loops"
 ```
 for \*nix shells,
 ```
 $env:FPM_FFLAGS="-ffree-line-length-none -fPIC -w -fallow-argument-mismatch -O3 -funroll-loops"
 ``` 
 for the PowerShell and
 ```
 set FPM_FFLAGS="-ffree-line-length-none -fPIC -w -fallow-argument-mismatch -O3 -funroll-loops"
 ```
 for the windows cmd prompt.

Authors and contributors  
========================

+   [Samuele Giuli](https://github.com/SamueleGiuli)  
    +   PhD student in Condensed Matter Physics @ SISSA, Trieste, Italy


+   [Gabriele Bellomia](https://github.com/beddalumia)  
    +   Post-Doc in Condensed Matter Physics @ SISSA, Trieste, Italy


+   [Carlos Mejuto Zaera](https://github.com/CarlosMejZ)  
    +   Post-Doc in Condensed Matter Physics @ SISSA, Trieste, Italy

	
