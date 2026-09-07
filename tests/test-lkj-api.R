suppressPackageStartupMessages(library(INLAcircular))

# LKJ.eta has no effect and is not validated when no special term is present.
ordinary <- LKJcc(y ~ 1, LKJ.eta = 0)
stopifnot(
  identical(ordinary$lkj_vars, character()),
  identical(ordinary$models, list())
)

if (requireNamespace("graphpcor", quietly = TRUE) &&
    requireNamespace("INLAtools", quietly = TRUE)) {
  prepared <- suppressWarnings(
    LKJcc(
      list(
        y1 ~ f(i, model = "iidkd_LKJ"),
        y2 ~ f(i, model = "iidkd_LKJ")
      ),
      LKJ.eta = 2,
      n.obs = 5
    )
  )

  process <- prepared$processes$i
  stopifnot(
    identical(process$pc.prior.u, c(1, 1)),
    identical(process$pc.prior.alpha, c(0.5, 0.5)),
    length(prepared$models) == 1L,
    length(prepared$formulas) == 2L
  )

  formula_text <- paste(
    vapply(
      prepared$formulas,
      function(formula) paste(deparse(formula), collapse = " "),
      ""
    ),
    collapse = " "
  )
  stopifnot(
    !grepl("pc.prior.u", formula_text, fixed = TRUE),
    !grepl("pc.prior.alpha", formula_text, fixed = TRUE)
  )
}
