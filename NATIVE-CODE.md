# Native-code organization

The R-package build and the INLA source-tree build are deliberately separate.

```text
src/
  INLAcirc_distributions.c       R-facing VM, LAvM, and PC-prior functions
  INLAcirc_r_bessel.c            R-facing Bessel wrapper
  INLAcirc_init.c                R native-routine registration
  INLAcirc_cloglike.c            package-DLL copy of the LAvM cloglike
  INLAcirc_cloglike.h            package-DLL likelihood prototype
  INLAcirc_cgeneric_compat.h     package-local minimal cloglike ABI
  INLAcirc_common.h              package-local numerical functions

src-cloglike/
  INLAcirc_cloglike.c            standalone INLA likelihood
  INLAcirc_cloglike.h            prototype using INLA's cgeneric.h
  INLAcirc_common.h              standalone ISO C numerical functions
  tests/INLAcirc_cloglike_test.c native C regression test
  tests/test-cloglike-inla.R     R/INLA integration test
```

## Hard build boundary

`R CMD INSTALL` compiles translation units found in `src/` only. No file in
`src/` includes a source or header from `src-cloglike/`, and no file in
`src-cloglike/` is pulled into the package DLL. In particular, package
installation does not need INLA's `cgeneric.h`.

The standalone `src-cloglike/INLAcirc_cloglike.h` includes `cgeneric.h` by
default. That header is intentionally not distributed: it must be supplied by
the INLA source tree when the module is compiled into inla-build. The test-only
ABI declarations under `src-cloglike/tests/` are enabled solely for the local
regression tests.

The implementations of `INLAcirc_cloglike.c` and `INLAcirc_common.h` in the
two source trees are kept byte-for-byte identical. When either implementation
changes, copy the same revision to both locations. This small duplication is
intentional: it prevents either build system from crossing the runtime
boundary.

## Direct `inla()` integration

Loading `INLAcircular` after `INLA` exposes `INLAcircular::inla()` on the
search path. It delegates ordinary likelihoods unchanged to `INLA::inla()`.
For `family = "lavm"`, it translates the call to the package's compiled
`cloglike`, converts the response with `INLA::inla.mdata()`, and passes the
LAvM link, concentration, and PC-prior controls to the native module. This is
also the path used internally by `inlacc()`.

The public control is `hyper$kappa`. Its `initial` value is already on INLA's
internal `log(kappa)` scale; there is no separate `log.initial` argument.
Supported prior names are `pc.vm0` and `pc.vminf`.

## Namespace policy

All INLAcircular implementation symbols and helpers begin with `INLAcirc_`.
All preprocessor identifiers begin with `INLACIRC_`. The `inla_` type names
in the two compatibility headers are fixed names from INLA's external cloglike
ABI, not package-owned implementation symbols.

The exported likelihood symbol is:

```text
INLAcirc_cloglike_lavm
```

The R `.Call` symbols are similarly namespaced, for example
`INLAcirc_C_dvm`, `INLAcirc_C_dpc_vm0`, `INLAcirc_C_dpc_vminf`, and
`INLAcirc_C_bessel_i`.

## Numerical-interface correction

The split makes the Bessel normalizer's scale explicit:

- `INLAcirc_log_bessel_i0_scaled(kappa)` accepts natural-scale concentration.
- `INLAcirc_log_bessel_i0_scaled_from_log_kappa(log_kappa)` accepts the
  internal cloglike hyperparameter.

The R `dvm()` and `dlavm()` paths use the natural-scale function. The cloglike
path uses the log-scale function, so its intended parameterization is
unchanged.

## Verification

From the package root, run the R-free native regression test:

```sh
make -C src-cloglike check
make -C src-cloglike clean
```

With R and INLA installed, run the end-to-end standalone-source test:

```sh
Rscript src-cloglike/tests/test-cloglike-inla.R
```

The R test compiles `src-cloglike/INLAcirc_cloglike.c` into a temporary shared
library and passes it to `INLA::inla.cloglike.define()`. It therefore tests
the standalone source directly rather than the package DLL.
