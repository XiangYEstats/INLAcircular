# `graphpcor` LKJ integration

The special term

```r
f(i, model = "iidkd_LKJ")
```

declares one component of an LKJ-correlated random vector. `INLAcircular`
groups all `iidkd_LKJ` terms with the same first argument (`i` above), in
likelihood order. For a group occurring in `m` likelihoods, `LKJcc()` builds
the equivalent of

```r
cov_model <- cgeneric(
  "LKJ",
  n = m,
  eta = LKJ.eta,
  sigma.prior.reference = collected_pc_prior_u,
  sigma.prior.probability = collected_pc_prior_alpha
)
```

and rewrites each special term to an ordinary INLA term using that cgeneric
model and an automatically generated replicate index.

## User interface

- The shared LKJ shape `LKJ.eta` defaults to `5`. Override it once with
  `inlacc(..., LKJ.eta = value)` when needed.
- `pc.prior.u` and `pc.prior.alpha` default to `1` and `0.5`, respectively,
  in every special `f()` term. Supply either argument only when a component
  needs a different PC prior.
- Do not create `i` or a replicate column for the usual case. The component
  index is generated in LKJ-group order and the replicate index is
  `seq_len(nrow(data))`.
- To provide the replicate mapping explicitly, pass

  ```r
  latent.index = list(
    index(
      var = "i",
      data.id = your_replicate_ids,
      process.id = "likelihood"
    )
  )
  ```

  Here `data.id` is deliberately used for INLA's `replicate` argument;
  `process.id = "likelihood"` requests automatic component indices. This
  special `index()` object is consumed by `LKJcc()` and is not passed through
  the ordinary latent-index path.

## Implementation boundary

`R/cc.LKJ.R` owns detection, prior collection, validation, graphpcor model
construction, formula rewriting, and LKJ index creation. The LKJ-specific
work in `R/inlacc.R` is limited to calling `LKJcc()`, excluding its internal
columns from generic index generation, merging prepared block columns after
ordinary `f()` data, exposing the model objects in the formula environment,
and recording fitted LKJ modes/metadata.

The merge order is important. In the earlier draft, the generic `f()` data
copy overwrote the LKJ component index with `1:n`, so the covariance model no
longer received component values `1:m`.

## Matching a direct INLA fit

For a direct comparison, the likelihood hyperparameters and requested output
must also match. In particular,

```r
family.setting = list(
  link = "inverse.tangent",
  hyper = list(kappa = list(initial = 15, fixed = TRUE))
)
```

fixes the LAvM log-concentration at 15. A fixed LAvM likelihood contributes no
free hyperparameter; the six-dimensional LKJ component therefore contributes
all 21 free hyperparameters in the biomechanical comparison. Use
`metrics = FALSE` with `control.compute = list(config = TRUE)` to match the
standalone script's computation request. Setting `metrics = TRUE` additionally
computes CPO, DIC, and WAIC and is not a like-for-like timing comparison.

## Biomechanical example

The biomechanical regression test and its six CSV files are stored in the
source package under `tests/manual`. Open
`tests/manual/07-biomechanical-lkj.R` in RStudio and click **Source**, or run:

```sh
Rscript tests/manual/07-biomechanical-lkj.R
```

No working-directory setup is required: the script locates the data relative
to its own file. It reproduces six likelihood blocks: three LAvM responses,
three Gaussian responses, one six-dimensional LKJ component named `i`, and
one replicate per biomechanical observation.
