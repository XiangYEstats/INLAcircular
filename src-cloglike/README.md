# R-independent INLA cloglike source

This directory is the canonical implementation of the LAvM `cloglike` module.
It is deliberately independent of the R C API and can be compiled without
linking against R or Rmath.

Files:

- `INLAcirc_cloglike.c`: likelihood, PC-prior, and INLA command dispatcher.
- `INLAcirc_cloglike.h`: public prototype for `INLAcirc_cloglike_lavm()`.
- `INLAcirc_common.h`: shared Bessel and PC-prior numerical functions; ISO C
  and `libm` only.
- `INLAcirc_cgeneric_compat.h`: minimal standalone declaration of the INLA ABI.
- `tests/INLAcirc_cloglike_test.c`: R-free native regression tests.

Run the standalone check with:

```sh
make check
make clean
```

## Adding the module to inla-build

Compile `INLAcirc_cloglike.c` as an INLA translation unit and register the
symbol:

```text
INLAcirc_cloglike_lavm
```

Use the INLA source tree's own `cgeneric.h` by adding
`-DINLACIRC_USE_INLA_CGENERIC` and making that header available on the include
path. The compatibility header is only for standalone/package compilation.

Do not copy or compile `../src/INLAcirc_cloglike_bridge.c` in inla-build. That
bridge exists solely because the current R package needs the symbol in its own
shared library until the likelihood is available in the INLA executable.
