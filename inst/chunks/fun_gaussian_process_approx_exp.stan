  /* Spectral density function of a Gaussian process
   * with exponential covariance kernel
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns:
   *   numeric values of the function evaluated at 'x'
   */
  vector spd_cov_exp_quad(data array[] vector x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real constant = square(sdgp) * pow(2, D) * pow(pi(), (D - 1) / 2.0) * tgamma((D + 1) / 2.0) / lscale[1];
      real inverse_l2 = 1 / square(lscale[1]);
      real four_pi2 = 4 * square(pi());
      real power = -(D + 1) / 2.0;
      for (m in 1:NB) {
        out[m] = constant * pow(inverse_l2 + four_pi2 * dot_self(x[m]), power);
      }
    } else {
      // multi-dimensional non-isotropic GP
      real norm_l = sqrt(dot_self(lscale));
      real constant = square(sdgp) * pow(2, D) * pow(pi(), (D - 1) / 2.0) * tgamma((D + 1) / 2.0) / norm_l;
      real inverse_l2 = 1 / square(norm_l);
      real four_pi2 = 4 * square(pi());
      real power = -(D + 1) / 2.0;
      for (m in 1:NB) {
        out[m] = constant * pow(inverse_l2 + four_pi2 * dot_self(x[m]), power);
      }
    }
    return out;
  }
