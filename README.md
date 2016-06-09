# 2016-SIOPT-stratification

Matlab code to reproduce the results of the paper:

> Vincent Duval, Jalal Fadili, Jérôme Malick, Gabriel Peyré,
> Sensitivity Analysis for Mirror-Stratified Convex Functions,
> Preprint arxiv, 2016.

The main functions are, for * in {'l1' 'nuclear'}:
- test_*.m: main scripts to reproduce the figures (phase transition and FB paths).
- perform_*_reg_fb.m: solve the regularization using Forward-Backward (FB).
- perform_*_reg_cvx.m: solve the regularization using interior point (with CVX)
- compute_certificate_*.m: compute the minimal norm certificate (the one that govern support/rank stability).

_IMPORTANT:_ this code requires [CVX](cvxr.com) to be installed.

Copyright (c) 2016 Gabriel Peyre
