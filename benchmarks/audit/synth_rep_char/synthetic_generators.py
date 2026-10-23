"""SYNTHETIC generators for the replica-comparison characterization (REPCHAR-1). Not production code.

All series have exact, stated stationary marginals; temporal correlation never changes the marginal.

* ``ar1_latent(rng, n, phi)``: Z_0 ~ N(0, 1) (stationary start), Z_t = phi Z_{t-1} + sqrt(1 - phi^2) e_t, e_t ~ N(0, 1).
  Marginal N(0, 1) for every t; lag-k autocorrelation phi^k.
* Linear channel: X_t = mu + sigma Z_t (marginal N(mu, sigma^2)).
* Circular channel: theta_t = wrap(mu + s Z_t) to (-pi, pi] (marginal wrapped normal; circular SD = s exactly,
  mean resultant length exp(-s^2 / 2)).
* Two-state mixture: hidden Markov chain S_t in {0, 1} with P(0 -> 1) = p01 = w / D and P(1 -> 0) = p10 = (1 - w) / D,
  so the stationary weight of state 1 is w and the state autocorrelation is (1 - 1/D)^k. S_0 is drawn from the
  stationary distribution. Mean dwell times are D / w in state 0 and D / (1 - w) in state 1 (D = 1: independent states).
  Requires D >= 1 and 0 < w < 1 so that both probabilities are <= 1. Emission: X_t = m_{S_t} + sigma_c e_t, e_t iid N(0, 1).
"""

import numpy as np

TWO_PI = 2.0 * np.pi


def ar1_latent(rng, n, phi):
    if not (0.0 <= phi < 1.0):
        raise ValueError("phi must be in [0, 1)")
    e = rng.standard_normal(n)
    z = np.empty(n)
    z[0] = e[0]
    c = np.sqrt(1.0 - phi * phi)
    for t in range(1, n):
        z[t] = phi * z[t - 1] + c * e[t]
    return z


def linear_series(rng, n, mu, sigma, phi):
    return mu + sigma * ar1_latent(rng, n, phi)


def wrap(theta):
    return np.pi - np.mod(np.pi - np.asarray(theta, dtype=np.float64), TWO_PI)  # (-pi, pi]


def circular_series(rng, n, mu, s, phi):
    return wrap(mu + s * ar1_latent(rng, n, phi))


def markov_states(rng, n, w, D):
    if not (0.0 < w < 1.0) or D < 1.0:
        raise ValueError("need 0 < w < 1 and D >= 1")
    p01, p10 = w / D, (1.0 - w) / D
    u = rng.random(n)
    s = np.empty(n, dtype=np.int8)
    s[0] = 1 if u[0] < w else 0
    for t in range(1, n):
        s[t] = (1 if u[t] < p01 else 0) if s[t - 1] == 0 else (0 if u[t] < p10 else 1)
    return s


def mixture_series(rng, n, m0, m1, sigma_c, w, D):
    s = markov_states(rng, n, w, D)
    return np.where(s == 1, m1, m0) + sigma_c * rng.standard_normal(n), s


def ess_gaussian_mean(n, phi):
    """Gaussian sample-mean effective sample size n (1 - phi) / (1 + phi). An approximation for the mean of a Gaussian
    AR(1) series only; not a guarantee for histograms, JS distances or circular features."""
    return n * (1.0 - phi) / (1.0 + phi)
