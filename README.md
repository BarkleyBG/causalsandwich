
# `causalgeex` for estimating causal effects with empirical sandwich variance estimator

Currently implemented as of version 0.0.0.9000:

- Estimators of the average treament effect
- IPTW
- G-formula
- Doubly Robust estimator

## Why use `causalgeex`? Easy implementation with valid inference!

- This package is inteded to be a "one-stop-shop" with emphasis on ease of use.
- The estimators here are consistent and asymptotically normal following usual regularity assumptions.


### What's `geex`?

[geex](https://cran.r-project.org/web/packages/geex/index.html) is a package designed for easy implementation of estimating equations. The `causalgeex` package is powered by `geex`

### What are estimating equations?

See [Stefanski and Boos (2002)](http://www.jstor.org/stable/3087324?seq=1#page_scan_tab_contents)

## How to install

The package can be installed from Github:

`devtools::install_github("BarkleyBG/causalgeex")`
