## MIDACO wrappers

Contains basic, naive implementations of wrappers between the C version of MIDACO (supplied as an object file) and Fortran subroutines,
via an intermediary C function (`midaco.c`) to transfer license-keys.

To use in Fortran:

```Fortran

include midaco_interface

call midaco(p, o, n, ni, m, me, x, f, g, xl, xu, iflag, &
            istop, param, rw, lrw, iw, liw, pf, lpf, save2file, &
            maxeval, maxtime, printeval)
            
```

where the variables have the same data-type as in the Fortran MIDACO codes; this will handle both the reverse communicator call to
the optimiser, and also the call to the printing routines. At the moment, all of the Fortran `REAL`/`INTEGER`-type
variables are converted to `C`-type using the `ISO_C_BINDING` module on every function call. In the future, this will be reduced to only
`iflag` and `pf`, so that conversion overheads can be limited. Issues with using `TARGET` also prevent many compiler optimisations from being used.

An additional `licenceKey.h` file is required in the root directory, and must contain a declaration of the license key to be used.

```C

char key[] = "*******************************************************"

```

this key is `#include`d in the `midaco.c` C source code. `midaco.c` currently supports a naive implementation for determining whether
the MIDACO subroutine has been called before; this is to allow optimiser re-starts in parallel routines.
