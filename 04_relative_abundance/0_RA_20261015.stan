data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1, upper=T> survey;
  array[N] real<lower=0, upper=1> y;
}

parameters {
  vector[T] alpha;
  real log_phi;
}

transformed parameters {
  real<lower=0> phi;
  vector<lower=0, upper=1>[T] mu_raw;
  vector<lower=0, upper=1>[T] mu;

  phi = exp(log_phi) + 1e-6;
  mu_raw = inv_logit(alpha);

  for (t in 1:T) {
    mu[t] = fmin(fmax(mu_raw[t], 1e-9), 1 - 1e-9);
  }
}

model {
  alpha ~ normal(0, 1.5);
  log_phi ~ normal(0, 1);

  for (n in 1:N) {
    y[n] ~ beta(mu[survey[n]] * phi,
                (1 - mu[survey[n]]) * phi);
  }
}

generated quantities {
  vector<lower=0, upper=1>[T] mu_out;
  mu_out = mu;
}
