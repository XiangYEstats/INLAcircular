## Test environments

- Fedora Linux 44, x86_64, R 4.6.1 (R-release), GCC 16.2.1:
  0 errors, 0 warnings, 3 notes.
- Fedora Linux 44 R-hub gcc16 container, x86_64, R-devel r90473, GCC 16.2.1
  (`--no-manual`): 0 errors, 0 warnings, 1 note.
- Ubuntu 20.04.6, x86_64, R 4.1.3, hard dependencies only and no INLA:
  installation and all regular tests passed.
- Standalone native tests under strict ISO C17 with GCC 16.2.1 and Clang
  22.1.8: passed.
- Full package source installation and compiled native-call smoke test with
  Clang 22.1.8 on Fedora: passed.
- Current stable INLA 26.08.07, including all seven developer model scripts
  under `tests/manual`: passed.

Windows, macOS ARM64/Intel, and Ubuntu R-release/R-devel are configured in
GitHub Actions and will be reported here after those remote checks have
completed.

## R CMD check results

There were no errors or warnings. There were three notes:

1. This is a new submission. INLA is in Suggests and is available from the
   repository declared in `Additional_repositories`.
2. The Fedora R-devel installation injected compiler hardening and
   architecture flags. The package has no Makevars file and does not add those
   flags itself.
3. HTML manual math-rendering checks were skipped because V8 was unavailable
   in the local check environment. The package HTML, vignette rebuild, and PDF
   manual checks otherwise completed successfully.

## Optional dependency behavior

INLA is intentionally suggested rather than imported so that the core
distribution and prior functions install without it. Model-fitting entry
points check for INLA >= 25.08.21 and the required exported cloglike API, then
provide an actionable installation or upgrade message. The hard-dependencies-
only, no-INLA test profile passes.
