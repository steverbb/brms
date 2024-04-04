  /* Spectral density function of a Gaussian process
   * with exponential covariance kernel
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns:
   *   numeric values of the function evaluated at 'x'
   */
  vector spd_cov_exp(data array[] vector x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real constant = square(sdgp) * lscale[1])^D; # edit it
      real neg_half_lscale2 = -0.5 * square(lscale[1]); # edit it
      for (m in 1:NB) {
        out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m])); # edit it
      }
    } else {
      // multi-dimensional non-isotropic GP
      real constant = square(sdgp) * prod(lscale); # edit it
      vector[Dls] neg_half_lscale2 = -0.5 * square(lscale); # edit it
      for (m in 1:NB) {
        out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m]))); # edit it
      }
    }
    return out;
  }
