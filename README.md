# libjupitermag

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.7306035.svg)](https://zenodo.org/badge/latestdoi/560414844)

Code for obtaining magnetic field vectors and traces from within Jupiter's magnetosphere using various magnetic field models.

This is part of a community code project :

[Magnetospheres of the Outer Planets Group Community Code](https://lasp.colorado.edu/mop/missions/juno/community-code/)

**Journal Paper DOI**: [https://doi.org/10.1007/s11214-023-00961-3](https://doi.org/10.1007/s11214-023-00961-3)
(PDF via DOI, or [https://rdcu.be/c5I71](https://rdcu.be/c5I71), see [Journal Publication](README.md#journal-publication).)

**Authors**

- Matt James - University of Leicester

- Gabby Provan - University of Leicester

- Aneesah Kamran - University of Leicester

- Rob Wilson - LASP

- Marissa Vogt - Boston University

- Marty Brennan - NASA JPL

- Stan Cowley - University of Leicester

This module forms part of the [JupiterMag](https://github.com/mattkjames7/JupiterMag.git) package for Python.

## Cloning and Building

This module requires a few submodules to be fetched, so the following command should clone everything:

```bash
git clone --recurse-sumodules https://github.com/mattkjames7/libjupitermag.git
```

This library requires `g++`, `make` and `ld` (Linux) or `libtool` (Mac) in order to be compiled. On Windows these tools can be provided by TDM-GCC.

To build in Linux and Mac simply run

```bash
cd libjupitermag
make
#optionally install the library
sudo make install
```

where the installation defaults to `/usr/local` but can be changed using the `PREFIX` argument, e.g.:

```bash
sudo make install PREFIX=/usr
```

It can also be uninstalled,

```bash
make uninstall
```

Under Windows powershell/command line,

```powershell
cd libjupitermag
./compile.bat
```

or under Linux but building for Windows

```bash
cd libjupitermag
make windows
```

After a successful build, a new library (`libjupitermag.so`) or DLL (`libjupitermag.dll`) should appear in the `lib/libjupitermag` directory.

## Usage

The shared object library created by compiling this project should be accessible by various programming languages. This section mostly covers C++, other languages should be able to use the functions defined in the `extern "C"` section of the header file quite easily. For Python, it is straightforward to use `ctypes` to access this library, similarly IDL could use the `CALL_EXTERNAL` function.

### Linking to the library

Here is a very basic example of how to link to the code using C++ and print the version of the library:

```cpp
/* test.cc */
#include <stdio.h>
#include <jupitermag.h>

int main() {
    /* simply print the version of the library */
    printf("libjupitermag version: $%d.%d.%d\n",
            LIBJUPITERMAG_VERSION_MAJOR,
            LIBJUPITERMAG_VERSION_MINOR,
            LIBJUPITERMAG_VERSION_PATCH);
    return 0;
} 
```

If the library was installed using `sudo make install` then the library can be linked to during compiling time with:

```bash
g++ test.c -o test -ljupitermag
```

Otherwise, the header should be included using a relative or absolute path, e.g:

```cpp
/* instead of this */
#include <jupitermag>
/* use seomthing like this */
#include "include/jupitermag.h"
```

then compile, e.g.

```gcc
#absoulte path
g++ test.cc -o test -L:/path/to/libjupitermag.so

#or relative path
g++ test.cc -o test -Wl,-rpath='$$ORIGIN/../lib/libjupitermag' -L ..lib/libjupitermag -ljupitermag
```

### Calling Field Models

Internal field models can be called using the `InternalModel` object, e.g.

```cpp
#include <stdio.h>
#include <jupitermag.h>

int main () {
    /* create an instance of the object */
    InternalModel modelobj = InternalModel();

    /* set the model to use */
    modelobj.SetModel("jrm09");

    /* get the model vectors at some position */
    double x = 10.0, y = 0.0, z = 0.0;
    double Bx, By, Bz;
    modelobj.Field(x,y,z,&Bx,&By,&Bz);

    printf("B at [%f,%f,%f] = [%f,%f,%f]\n",x,y,z,Bx,By,Bz);
}
```

There are also simple functions for each of the models included in the library, see the table below.

| Model              | C String    | Field Function   | Reference                                                            |
| ------------------ | ----------- | ---------------- | -------------------------------------------------------------------- |
| JRM33              | `jrm33`     | `jrm33Field`     | Connerney et al., 2022                                               |
| JRM09              | `jrm09`     | `jrm09Field`     | Connerney et al., 2018                                               |
| ISaAC              | `isaac`     | `isaacField`     | Hess et al., 2017                                                    |
| VIPAL              | `vipal`     | `vipalField`     | Hess et al., 2011                                                    |
| VIP4               | `vip4`      | `vip4Field`      | Connerney 2007                                                       |
| VIT4               | `vit4`      | `vit4Field`      | Connerney 2007                                                       |
| O4                 | `o4`        | `o4Field`        | Connerney 1981                                                       |
| O6                 | `o6`        | `o6Field`        | Connerney 2007                                                       |
| GSFC15evs          | `gsfc15evs` | `gsfc15evsField` | Connerney 1981                                                       |
| GSFC15ev           | `gsfc15ev`  | `gsfc15evField`  | Connerney 1981                                                       |
| GSFC13ev           | `gsfc13ev`  | `gsfc13evField`  | Connerney 1981                                                       |
| Ulysses 17ev       | `u17ev`     | `u17evField`     | Connerney 2007                                                       |
| SHA                | `sha`       | `shaField`       | Connerney 2007                                                       |
| Voyager 1 17ev     | `v117ev`    | `v117evField`    | Connerney 2007                                                       |
| JPL15ev            | `jpl15ev`   | `jpl15evField`   | Connerney 1981                                                       |
| JPL15evs           | `jpl15evs`  | `jpl15evsField`  | Connerney 1981                                                       |
| P11A               | `p11a`      | `p11aField`      |                                                                      |
| **External Model** |             |                  |                                                                      |
| Con 2020           | `con2020`   | `Con2020Field`   | Connerney et al., 1981; Edwards et al., 2001; Connerney et al., 2020 |

### Field Tracing

Field tracing can be done using the `Trace` object, e.g.

```cpp
#include <stdio.h>
#include <jupitermag.h>
#include <vector>

int main () {
    /* set initial position to start trace from (this can be an array
        for multiple traces) */
    int n = 1;
    double x0 = 5.0;
    double y0 = 0.0;
    double z0 = 0.0;
    int nalpha = 1;
    double alpha = 0.0;

    printf("Create field function vector\n");
    /* store the function pointers for the components of the
    model to be included in the trace */
    std::vector<FieldFuncPtr> Funcs;

    /* internal model */
    Funcs.push_back(jrm09Field);

    /* external model */
    Funcs.push_back(Con2020Field);

    /* initialise the trace object */
    printf("Create Trace object\n");
    Trace T(Funcs);

    /* add the starting posiutions fo the traces */
    printf("Add starting position\n");
    T.InputPos(n,&x0,&y0,&z0);

    /* configure the trace parameters, leaving this empty will
    use default values for things like minimum and maximum step size  */
    printf("Set the trace parameters \n");
    T.SetTraceCFG();

    /* set up the alpha calculation - the angles (in degrees) of each 
    polarization angle. This is generally used for ULF waves */
    printf("Initialize alpha\n");
    T.SetAlpha(nalpha,&alpha);    

    /* Trace */
    printf("Trace\n");
    T.TraceField();

    /* trace distance, footprints, Rnorm */
    printf("Footprints etc...\n");
    T.CalculateTraceDist();
    T.CalculateTraceFP();
    T.CalculateTraceRnorm();

    /* calculate halpha for each of the polarization angles 
       specified above*/
    printf("H_alpha\n");
    T.CalculateHalpha();
}
```

The above code traces along the magnetic field using the JRM09 internal and Con2020 external field models together. The trace coordinates and field vectors at each step can be obtained from the member variables `T.x_`, `T.y_`, `T.z_`, `T.Bx_`, `T.By_` and `T.Bz_`, where each is a 2D array with the shape (`T.n_`,`T.MaxLen_`), where `T.n_` is the number of traces and `T.MaxLen_` is the maximum number of steps allowed in the trace. The number of steps takin in each trace is defined in the `T.nstep_` array.

## References

- Connerney, J. E. P., Timmins, S., Herceg, M., & Joergensen, J. L. (2020). A Jovian magnetodisc model for the Juno era. *Journal of Geophysical Research: Space Physics*, 125, e2020JA028138. https://doi.org/10.1029/2020JA028138

- Connerney, J. E. P., Acuña, M. H., and Ness, N. F. (1981), Modeling the Jovian current sheet and inner magnetosphere, *J. Geophys. Res.*, 86( A10), 8370– 8384, doi:[10.1029/JA086iA10p08370](https://doi.org/10.1029/JA086iA10p08370).

- Connerney, J. E. P. (1981), The magnetic field of Jupiter: A generalized inverse approach, *J. Geophys. Res.*, 86( A9), 7679– 7693, doi:[10.1029/JA086iA09p07679](https://doi.org/10.1029/JA086iA09p07679 "Link to external resource: 10.1029/JA086iA09p07679").

- Connerney, J. E. P., Acuña, M. H., and Ness, N. F. (1982), Voyager 1 assessment of Jupiter's planetary magnetic field, *J. Geophys. Res.*, 87( A5), 3623– 3627, doi:[10.1029/JA087iA05p03623](https://doi.org/10.1029/JA087iA05p03623 "Link to external resource: 10.1029/JA087iA05p03623").

- Connerney, J.E.P.. (2007). Planetary Magnetism. Treatise on Geophysics. 10. 243-280. 10.1016/B978-044452748-6.00159-0.

- Connerney, J. E. P., Kotsiaros, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Joergensen, P. S., et al. (2018). A new model of Jupiter's magnetic field from Juno's first nine orbits. Geophysical Research Letters, 45, 2590– 2596. [A New Model of Jupiter's Magnetic Field From Juno's First Nine Orbits - Connerney - 2018 - Geophysical Research Letters - Wiley Online Library](https://doi.org/10.1002/2018GL077312)

- Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's Prime Mission. Journal of Geophysical Research: Planets, 127, e2021JE007055. [A New Model of Jupiter's Magnetic Field at the Completion of Juno's Prime Mission - Connerney - 2022 - Journal of Geophysical Research: Planets - Wiley Online Library](https://doi.org/10.1029/2021JE007055)

- Hess, S. L. G., Bonfond, B., Zarka, P., and Grodent, D. (2011), Model of the Jovian magnetic field topology constrained by the Io auroral emissions, *J. Geophys. Res.*, 116, A05217, doi:[10.1029/2010JA016262](https://doi.org/10.1029/2010JA016262 "Link to external resource: 10.1029/2010JA016262").Redirecting

- Hess, S., Bonfond, B., Bagenal, F., & Lamy, L. (2017). A model of the Jovian internal field derived from in-situ and auroral constraints, doi:[10.1553/PRE8s157](https://doi.org/10.1553/PRE8s157)

- Edwards T.M., Bunce E.J., Cowley S.W.H. (2001), A note on the vector potential of Connerney et al.'s model of the equatorial current sheet in Jupiter's magnetosphere, *Planetary and Space Science,*49, 1115-1123,[https://doi.org/10.1016/S0032-0633(00)00164-1](https://doi.org/10.1016/S0032-0633(00)00164-1).
