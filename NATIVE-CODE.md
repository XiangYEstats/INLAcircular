# Native-code organization

The native sources are separated by their host runtime.

```text
src/
  INLAcirc_distributions.c       R-facing distribution functions
  INLAcirc_r_bessel.c            R-facing Bessel wrapper
  INLAcirc_init.c                R native-routine registration
  INLAcirc_cloglike_bridge.c     package-only compatibility bridge

src-cloglike/
  INLAcirc_cloglike.c            R-independent INLA likelihood
  INLAcirc_cloglike.h            public likelihood prototype
  INLAcirc_common.h              shared ISO C numerical functions
  INLAcirc_cgeneric_compat.h     standalone INLA ABI declarations
```

## Dependency boundary

Files in `src/` may include R headers. No file in `src-cloglike/` may include
an R header or call the R/Rmath API. The reusable numerical code is defined as
`static inline` functions in `INLAcirc_common.h`, which includes only `math.h`.
Consequently, compiling the likelihood into inla-build requires `libm`, but
does not require R.

The compatibility bridge includes the canonical cloglike translation unit so
that existing calls to `INLA::inla.cloglike.define()` continue to work with
the package shared library. This dependency is one-way: the R package embeds
the pure module, while the pure module never depends on the R package.

## Namespace policy

All INLAcircular implementation symbols and helpers begin with `INLAcirc_`.
All preprocessor identifiers begin with `INLACIRC_`. The `inla_` type names in
`INLAcirc_cgeneric_compat.h` are fixed names from INLA's external cloglike ABI,
not package-owned implementation symbols.

The exported likelihood symbol is:

```text
INLAcirc_cloglike_lavm
```

The R `.Call` symbols are similarly namespaced, for example
`INLAcirc_C_dvm` and `INLAcirc_C_bessel_i`.

## Numerical-interface correction

The split makes the Bessel normalizer's scale explicit:

- `INLAcirc_log_bessel_i0_scaled(kappa)` accepts natural-scale concentration.
- `INLAcirc_log_bessel_i0_scaled_from_log_kappa(log_kappa)` accepts the
  internal cloglike hyperparameter.

Previously, the R density functions passed natural-scale `kappa` to a helper
that interpreted its argument as `log(kappa)`. The R `dvm()` and `dlavm()`
paths now use the natural-scale function. The cloglike path continues to use
the log-scale function, so its intended parameterization is unchanged.

## Verification

From the package root:

```sh
make -C src-cloglike check
make -C src-cloglike clean
```

The check compiles without R, exercises estimated and fixed concentration,
all three links, the prior command, and the CDF command, and confirms that the
object exposes no unprefixed implementation symbol.
