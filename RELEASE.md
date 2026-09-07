# Maintainer release workflow

`DESCRIPTION` must contain one static package version. Increase it once for
each CRAN submission that represents a new release; do not change it for every
commit.

Every update submitted after a version has been accepted on CRAN must have a
higher version. CRAN also prefers increasing the version for every
resubmission, even when the previous submission was not accepted, because it
reduces review and artifact confusion. Do not bump for ordinary development
commits that are never submitted.

Use semantic versioning for release numbers:

- patch: backward-compatible fixes, for example `1.0.0` to `1.0.1`;
- minor: backward-compatible features, for example `1.0.1` to `1.1.0`;
- major: intentional breaking changes, for example `1.4.0` to `2.0.0`.

During development after a release, a fourth component such as `1.0.0.9000`
may identify the unreleased branch. Do not submit a `.9000` development
version to CRAN.

## Prepare a release

### First-submission audit

Before the first CRAN submission:

1. Confirm the copyright and license provenance of every native-code fragment,
   especially the Acklam approximation, GSL-derived coefficient tables, INLA
   reference approximations, and cgeneric compatibility declarations. Preserve
   required upstream notices, add contributors/copyright holders to
   `Authors@R` or `inst/COPYRIGHTS` where appropriate, and use a compatible
   package license. This requires an authorship decision and must not be
   guessed from the code alone.
2. Test against the current stable INLA release from
   `https://inla.r-inla-download.org/R/stable`. Do not advertise the testing
   channel in a CRAN release. Keep INLA in `Suggests`, list its repository in
   `Additional_repositories`, require INLA >= 25.08.21, and guard all uses with
   the package's availability and compatibility checks.
3. Verify the C sources on Windows, macOS Intel/ARM, and Linux. In particular,
   retain the guarded fallbacks for `M_PI`, `M_LN2`, and `M_SQRT1_2`, and run
   the strict-C17 numerical assertions on every target.
4. Confirm that the graphpcor integration works with its current public
   release. The present integration accesses `cgeneric_LKJ` from the package
   namespace; prefer an exported API or obtain an explicit compatibility
   commitment from its maintainer before depending on that internal object.
5. Add an author-year citation with DOI/arXiv identifier to `Description` if
   this package implements a published method. Do not invent a citation.

### Release files

1. Choose the release type and update `Version:` in `DESCRIPTION`. The
   development helper `usethis::use_version("patch")` (or `"minor"` or
   `"major"`) can perform the version transition; editing `DESCRIPTION`
   manually is also valid.
2. Add a matching heading and user-visible changes to `NEWS.md`.
3. Regenerate Rd files with `roxygen2::roxygenise()` and confirm that
   `R CMD check` reports no code/documentation mismatch.
4. Re-knit `vignettes/INLAcircular-guide.Rmd` and inspect both repository
   outputs. The CRAN package vignette is HTML; the additional PDF under
   `output/pdf` is intentionally excluded from the source tarball.
5. Put release-specific check results and explanations in `cran-comments.md`.
   That file is maintainer correspondence and remains excluded by
   `.Rbuildignore`.

### Test the exact source tarball

From the package root, replace `X.Y.Z` with the release version:

```sh
R CMD build .
R CMD check --as-cran INLAcircular_X.Y.Z.tar.gz
```

- Build with the current R release or R-patched, and run the decisive check
  with current R-devel when possible.
- Do not use `--no-manual` for the release check.
- Require zero errors, zero warnings, and no unexplained significant notes.
- Run the native tests in `src-cloglike` separately without changing that
  directory's source files.
- Inspect `tar -tzf INLAcircular_X.Y.Z.tar.gz`; it must not contain compiled
  objects, `src-cloglike`, repository output, credentials, check directories,
  or other development files.
- Supplement the local check with Win-builder and a macOS check service (and
  R-hub when useful). This is essential for a package containing C code.
- Also test with suggested packages absent where practical; optional INLA and
  graphpcor code must skip cleanly.

### Automated checks in GitHub

The workflow in `.github/workflows/R-CMD-check.yaml` runs automatically for
pushes and pull requests. It checks R-release on Ubuntu, Windows, macOS ARM64,
and macOS Intel; R-devel on Ubuntu and Fedora; and a hard-dependencies-only
installation on the declared minimum R 4.1. It also compiles the standalone C
code in strict C17 mode and checks the current stable INLA integration.

Before each CRAN release, run the same workflow manually from **Actions > R CMD
check > Run workflow**. A manual run additionally executes all seven scripts
under `tests/manual` and uploads their logs. Do not submit until every required
job is green. GitHub checks supplement, but do not replace, the final local
`R CMD check --as-cran` on the exact source archive or the external CRAN
builder checks.

### Submit and follow up

1. Upload the exact checked `INLAcircular_X.Y.Z.tar.gz` source tarball through
   the CRAN web submission form; CRAN does not accept maintainer-built binary
   packages.
2. Accept the confirmation email. Do not upload another version while that
   submission is pending.
3. For a requested resubmission, explain each correction in the form's
   `Optional comment` field and normally increase the version.
4. After CRAN accepts the release, wait for the check page to finish updating,
   tag that version in version control, and optionally move the development
   branch to the corresponding `.9000` version.

For the next accepted release, repeat this entire workflow with a higher
semantic version. GitHub and CRAN are separate publication channels: push the
source and CI configuration to GitHub first, then upload the exact green source
archive to CRAN and complete its email confirmation.

Never compute the package version dynamically from the date, Git state, or
environment variables: CRAN metadata and installed package metadata must be
reproducible from `DESCRIPTION`.
