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

Authors and contributors  
========================

+   [Samuele Giuli](https://github.com/SamueleGiuli)  
    +   PhD student in Condensed Matter Physics @ SISSA, Trieste, Italy


+   [Gabriele Bellomia](https://github.com/beddalumia)  
    +   Post-Doc in Condensed Matter Physics @ SISSA, Trieste, Italy

	
