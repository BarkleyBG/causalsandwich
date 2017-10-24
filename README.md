[![Travis-CI Build Status](https://travis-ci.org/BarkleyBG/causalsandwich.svg?branch=master)](https://travis-ci.org/BarkleyBG/causalsandwich)
[![Coverage Status](https://img.shields.io/codecov/c/github/BarkleyBG/causalsandwich/master.svg)](https://codecov.io/github/BarkleyBG/causalsandwich?branch=master)

# `causalsandwich` is used to estimate causal effects with empirical sandwich variance estimator. This means that the estimators have nice asymptotic properties when models are correctly specified.

Currently implemented as of version 0.0.0.9000:

- Estimators of the average treament effect
- IPTW
- G-formula
- Doubly Robust estimator

## Why use `causalsandwich`? Easy implementation with valid inference!

- This package is inteded to be a "one-stop-shop" with emphasis on ease of use.
- The estimators here are consistent and asymptotically normal following usual regularity assumptions.


### What's `geex`?

[geex](https://cran.r-project.org/web/packages/geex/index.html) is a package designed for easy implementation of estimating equations. The `causalsandwich` package is powered by `geex`

### What are estimating equations?

See [Stefanski and Boos (2002)](http://www.jstor.org/stable/3087324?seq=1#page_scan_tab_contents)

## How to install

The package can be installed from Github:

`devtools::install_github("BarkleyBG/causalsandwich")`
