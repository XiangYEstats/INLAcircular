# Standalone INLA cloglike source

This directory contains the R-independent LAvM `cloglike` module that can be
compiled directly into the INLA source tree. Nothing here includes the R C API
or calls R/Rmath.

Files:

- `INLAcirc_cloglike.c`: likelihood, PC-prior, and INLA command dispatcher.
- `INLAcirc_cloglike.h`: public prototype for
  `INLAcirc_cloglike_lavm()`; it uses the INLA source tree's `cgeneric.h`.
- `INLAcirc_common.h`: Bessel and PC-prior numerical functions; ISO C and
  `libm` only.
- `FUNCTIONS`: exported function name used by the INLA build integration.
- `tests/INLAcirc_cloglike_test.c`: R-free native regression test.
- `tests/test-cloglike-inla.R`: end-to-end R test that calls INLA with this
  source compiled as a temporary shared library.

## Adding the module to inla-build

Place INLA's current `cgeneric.h` on the compiler include path (or add it to
this directory during the external build), compile `INLAcirc_cloglike.c`, and
register the symbol:

```text
INLAcirc_cloglike_lavm
```

`cgeneric.h` is intentionally absent from this package and is ignored by Git.
The ordinary production target therefore works only after the external header
has been supplied:

```sh
make all
```

The R package never compiles this directory. Its DLL is built solely from the
self-contained sources under `../src/`.

## Tests

The native regression test uses a test-only ABI declaration and does not need
R, INLA, or a copied `cgeneric.h`:

```sh
make check
make clean
```

With R and INLA installed, the integration test follows INLA's normal
`cloglike` shared-library workflow:

```sh
Rscript tests/test-cloglike-inla.R
```

It compiles `INLAcirc_cloglike.c` in a temporary directory, creates the
cloglike object with `INLA::inla.cloglike.define()`, and fits a small model
with `INLA::inla()`. The temporary build is removed automatically.
