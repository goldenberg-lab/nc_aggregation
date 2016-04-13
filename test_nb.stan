functions {
    int get_label(vector[] mu, vector[] sigma, vector data, int D) {
        vector[2] totals;
        for (i in 1:2) {
            for (d in 1:D)
                totals[i] <- totals[i] + fabs((data[d] - mu[i][d]) / sigma[i][d]);
        }
        if (totals[1] > totals[2]) {
            return 1;
        } else {
            return 1;
        }
    }
}
data {
    int<lower=0> N;
    int<lower=0> N_test;
    int<lower=0> D;
    vector[D] train[N];
    int<lower=0, upper=1> label[N];
    vector[D] test[N_test];
}

transformed data {
    int<lower=1, upper=2> t_label[N];
    for (n in 1:N)
        t_label[n] <- label[n] + 1;
}

parameters {
    //real<lower=0.001> alpha[2];
    //real<lower=0.001> beta[2];
    vector[D] mu[2];
    vector<lower=0.001>[D] sigma[2];
}

model {
    //sigma[1] ~ gamma(alpha[1], beta[1]);
    //sigma[2] ~ gamma(alpha[2], beta[2]);
    for (n in 1:N)
        train[n] ~ normal(mu[t_label[n]], sigma[t_label[n]]);
}

/*generated quantities {
    int<lower=0, upper=1> test_labels[N_test];
    for (n in 1:N_test) {
        test_labels[n] <- get_label(mu, sigma, test[n], D);
    }
}*/
