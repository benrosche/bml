# Changelog

## bml 0.9.0

### Major Changes

- **Package renamed from `rmm` to `bml`** (Bayesian Multiple-Membership
  Multilevel Models)

- **New syntax for weight functions:** The `ar` parameter has been moved
  from the [`fn()`](https://benrosche.github.io/bml/reference/fn.md)
  specification to the
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) block level
  for clearer API

  - Old: `fn(w ~ 1/n, c = TRUE, ar = FALSE)`
  - New: `fn(w ~ 1/n, c = TRUE)` with `ar = FALSE` at the
    [`mm()`](https://benrosche.github.io/bml/reference/mm.md) level

- **Support for multiple mmid groups:** The package now supports models
  with multiple membership identifiers, allowing more complex membership
  structures

- **Enhanced documentation:** Comprehensive documentation added for the
  `coalgov` dataset including:

  - Detailed variable descriptions organized by level (identifiers,
    government-level, country-level, party-level)
  - Statistical summaries for all variables
  - Clear explanation of multiple-membership structure
  - Updated references and examples

### New Features

- **Flexible weight function parameterization:** Enhanced support for
  parameterizing weight functions with covariates and group-specific
  structures

- **Per-group random effects:** Random effects can now be specified
  separately for different mmid groups

- **Improved JAGS code generation:** Optimized model string generation
  for better performance with complex multiple-membership structures

### Breaking Changes

- **`ar` parameter moved:** Existing code using `fn(w ~ ..., ar = TRUE)`
  must be updated to place `ar` in the
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) block
  instead

- **Dataset changes:**

  - Removed `schoolnets` dataset (including `nodedat` and `edgedat`
    objects)
  - Updated `coalgov` dataset with enhanced documentation and additional
    variables

### Bug Fixes

- Fixed issues with weight function constraints when using multiple
  [`mm()`](https://benrosche.github.io/bml/reference/mm.md) blocks

- Improved handling of group-level indices in JAGS variable creation

### Documentation

- Updated vignette examples to use new syntax
- Added comprehensive FAQ section
- Improved installation instructions with CRAN and GitHub options
- Updated references to reflect 2026 publication

------------------------------------------------------------------------
