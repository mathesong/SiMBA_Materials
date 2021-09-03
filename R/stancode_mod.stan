// generated with brms 2.14.4
functions {
  /* turn a vector into a matrix of defined dimension 
   * stores elements in row major order
   * Args: 
   *   X: a vector 
   *   N: first dimension of the desired matrix
   *   K: second dimension of the desired matrix 
   * Returns: 
   *   a matrix of dimension N x K 
   */ 
  matrix as_matrix(vector X, int N, int K) { 
    matrix[N, K] Y; 
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
    }
    return Y; 
  } 
 /* compute correlated group-level effects
  * Args: 
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns: 
  *   matrix of scaled group-level effects
  */ 
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }


real twotcm_log_stan(real logK1, real logVnd, real logBPnd, 
                    real logk4, real logvB, real time,
                    real t0_aif, real b_aif, real lambda1_aif, 
                    real lambda2_aif, real lambda3_aif, 
                    real A1_aif, real A2_aif, real A3_aif, 
                    real tstar_aif, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real Vnd;
    real BPnd;
    real vB;
    
    real R1;
    real R2;
    real L1;
    real L2;
    
    int indicatort0;
    int indicatorpeak;
    
    real timesincet0;
    real tstaraftert0;
    
    real pred;
    
    indicatort0 = time > t0_aif;
    indicatorpeak = time > tstar_aif;
    

    
    K1 = exp(logK1);
    Vnd = exp(logVnd);
    BPnd = exp(logBPnd);
    k4 = exp(logk4);
    vB = exp(logvB);
    
    k2 = K1 / Vnd;
    k3 = k4 * BPnd;
  
    R1 = 0.5 * (k2 + k3 + k4 + sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4));
    R2 = 0.5 * (k2 + k3 + k4 - sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4));

    L1 = (K1 * (R1 - k3 - k4))/(R1 - R2);
    L2 = (K1 * (k3 + k4 - R2))/(R1 - R2);
    
    timesincet0 = time - t0_aif;
    tstaraftert0 = tstar_aif - t0_aif;
  
  pred =  (indicatort0)*((1-vB)*  // whole thing is zero, if before t0
      
       ((-1*(indicatorpeak-1))*(b_aif*L1*(exp(-R1*timesincet0)/R1^2 + timesincet0/R1 - 1/R1^2) +
       b_aif*L2*(exp(-R2*timesincet0)/R2^2 + timesincet0/R2 - 1/R2^2)) +
       
       (indicatorpeak)*(b_aif*L1*exp(-R1*timesincet0)*(tstaraftert0/R1*exp(R1*tstaraftert0) - 
       1/R1^2*exp(R1*tstaraftert0) + 
       1/R1^2) + 
       L1*A1_aif*(exp(-lambda1_aif*timesincet0)/(R1-lambda1_aif) -
       exp(R1*tstaraftert0-R1*timesincet0-lambda1_aif*tstaraftert0)/(R1-lambda1_aif)) +
       L1*A2_aif*(exp(-lambda2_aif*timesincet0)/(R1-lambda2_aif) -
       exp(R1*tstaraftert0-R1*timesincet0-lambda2_aif*tstaraftert0)/(R1-lambda2_aif)) +
       L1*A3_aif*(exp(-lambda3_aif*timesincet0)/(R1-lambda3_aif) -
       exp(R1*tstaraftert0-R1*timesincet0-lambda3_aif*tstaraftert0)/(R1-lambda3_aif)) +
       b_aif*L2*exp(-R2*timesincet0)*(tstaraftert0/R2*exp(R2*tstaraftert0) - 
       1/R2^2*exp(R2*tstaraftert0) + 1/R2^2) +
       L2*A1_aif*(exp(-lambda1_aif*timesincet0)/(R2-lambda1_aif) -
       exp(R2*tstaraftert0-R2*timesincet0-lambda1_aif*tstaraftert0)/(R2-lambda1_aif)) +
       L2*A2_aif*(exp(-lambda2_aif*timesincet0)/(R2-lambda2_aif) -
       exp(R2*tstaraftert0-R2*timesincet0-lambda2_aif*tstaraftert0)/(R2-lambda2_aif)) +
       L2*A3_aif*(exp(-lambda3_aif*timesincet0)/(R2-lambda3_aif) -
       exp(R2*tstaraftert0-R2*timesincet0-lambda3_aif*tstaraftert0)/(R2-lambda3_aif)))) + 
       
       vB*bloodval);
       
  return(pred);
}

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_logK1;  // number of population-level effects
  matrix[N, K_logK1] X_logK1;  // population-level design matrix
  int<lower=1> K_logBPnd;  // number of population-level effects
  matrix[N, K_logBPnd] X_logBPnd;  // population-level design matrix
  int<lower=1> K_logVnd;  // number of population-level effects
  matrix[N, K_logVnd] X_logVnd;  // population-level design matrix
  int<lower=1> K_logk4;  // number of population-level effects
  matrix[N, K_logk4] X_logk4;  // population-level design matrix
  int<lower=1> K_logvB;  // number of population-level effects
  matrix[N, K_logvB] X_logvB;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  vector[N] C_4;
  vector[N] C_5;
  vector[N] C_6;
  vector[N] C_7;
  vector[N] C_8;
  vector[N] C_9;
  vector[N] C_10;
  vector[N] C_11;
  int<lower=1> K_sigma;  // number of population-level effects
  matrix[N, K_sigma] X_sigma;  // population-level design matrix
  // data for splines
  int Ks_sigma;  // number of linear effects
  matrix[N, Ks_sigma] Xs_sigma;  // design matrix for the linear effects
  // data for spline s(t_tac)
  int nb_sigma_1;  // number of bases
  int knots_sigma_1[nb_sigma_1];  // number of knots
  // basis function matrices
  matrix[N, knots_sigma_1[1]] Zs_sigma_1_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_sigma_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_sigma_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_logK1_1;
  vector[N] Z_3_logBPnd_2;
  vector[N] Z_3_logVnd_3;
  vector[N] Z_3_logk4_4;
  int<lower=1> NC_3;  // number of group-level correlations
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_logK1_1;
  vector[N] Z_4_logBPnd_2;
  vector[N] Z_4_logVnd_3;
  vector[N] Z_4_logk4_4;
  int<lower=1> NC_4;  // number of group-level correlations
  // data for group-level effects of ID 5
  int<lower=1> N_5;  // number of grouping levels
  int<lower=1> M_5;  // number of coefficients per level
  int<lower=1> J_5[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_5_logVnd_1;
  vector[N] Z_5_logk4_2;
  int<lower=1> NC_5;  // number of group-level correlations
  // data for group-level effects of ID 6
  int<lower=1> N_6;  // number of grouping levels
  int<lower=1> M_6;  // number of coefficients per level
  int<lower=1> J_6[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_6_logvB_1;
  // data for group-level effects of ID 7
  int<lower=1> N_7;  // number of grouping levels
  int<lower=1> M_7;  // number of coefficients per level
  int<lower=1> J_7[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_7_logvB_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_logK1] b_logK1;  // population-level effects
  vector[K_logBPnd] b_logBPnd;  // population-level effects
  vector[K_logVnd] b_logVnd;  // population-level effects
  vector[K_logk4] b_logk4;  // population-level effects
  vector[K_logvB] b_logvB;  // population-level effects
  vector[K_sigma] b_sigma;  // population-level effects
  vector[Ks_sigma] bs_sigma;  // spline coefficients
  // parameters for spline s(t_tac)
  // standarized spline coefficients
  // vector[knots_sigma_1[1]] zs_sigma_1_1;
  vector[knots_sigma_1[1]] rs_sigma_1_1;
  real<lower=0> sds_sigma_1_1;  // standard deviations of spline coefficients
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  // vector[N_1] z_1[M_1];  // standardized group-level effects
  vector[N_1] r_1;    // actual group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  // vector[N_2] z_2[M_2];  // standardized group-level effects
  vector[N_2] r_2;  // actual group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  // matrix[M_3, N_3] z_3;  // standardized group-level effects
  matrix[N_3, M_3] r_3;  // actual group-level effects
  cholesky_factor_corr[M_3] L_3;  // cholesky factor of correlation matrix
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  matrix[M_4, N_4] z_4;  // standardized group-level effects
  cholesky_factor_corr[M_4] L_4;  // cholesky factor of correlation matrix
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  // matrix[M_5, N_5] z_5;  // standardized group-level effects
  matrix[N_5, M_5] r_5;  // actual group-level effects
  cholesky_factor_corr[M_5] L_5;  // cholesky factor of correlation matrix
  vector<lower=0>[M_6] sd_6;  // group-level standard deviations
  vector[N_6] z_6[M_6];  // standardized group-level effects
  // vector[N_6] r_6;  // standardized group-level effects
  vector<lower=0>[M_7] sd_7;  // group-level standard deviations
  vector[N_7] r_7;  // actual group-level effects
}
transformed parameters {
  // actual spline coefficients
  vector[knots_sigma_1[1]] s_sigma_1_1;
  vector[N_1] r_1_sigma_1;  // actual group-level effects
  vector[N_2] r_2_sigma_1;  // actual group-level effects
  // matrix[N_3, M_3] r_3;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_logK1_1;
  vector[N_3] r_3_logBPnd_2;
  vector[N_3] r_3_logVnd_3;
  vector[N_3] r_3_logk4_4;
  matrix[N_4, M_4] r_4;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_4] r_4_logK1_1;
  vector[N_4] r_4_logBPnd_2;
  vector[N_4] r_4_logVnd_3;
  vector[N_4] r_4_logk4_4;
  // matrix[N_5, M_5] r_5;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_5] r_5_logVnd_1;
  vector[N_5] r_5_logk4_2;
  vector[N_6] r_6_logvB_1;  // actual group-level effects
  vector[N_7] r_7_logvB_1;  // actual group-level effects
  // compute actual spline coefficients
  s_sigma_1_1 = rs_sigma_1_1;
  // r_1_sigma_1 = (sd_1[1] * (z_1[1]));
  r_1_sigma_1 = r_1;
  // r_2_sigma_1 = (sd_2[1] * (z_2[1]));
  r_2_sigma_1 = r_2;
  // compute actual group-level effects
  // r_3 = scale_r_cor(z_3, sd_3, L_3);
  r_3_logK1_1 = r_3[, 1];
  r_3_logBPnd_2 = r_3[, 2];
  r_3_logVnd_3 = r_3[, 3];
  r_3_logk4_4 = r_3[, 4];
  // compute actual group-level effects
  r_4 = scale_r_cor(z_4, sd_4, L_4);
  r_4_logK1_1 = r_4[, 1];
  r_4_logBPnd_2 = r_4[, 2];
  r_4_logVnd_3 = r_4[, 3];
  r_4_logk4_4 = r_4[, 4];
  // compute actual group-level effects
  // r_5 = scale_r_cor(z_5, sd_5, L_5);
  r_5_logVnd_1 = r_5[, 1];
  r_5_logk4_2 = r_5[, 2];
  r_6_logvB_1 = (sd_6[1] * (z_6[1]));
  // r_6_logvB_1 = r_6;
  // r_7_logvB_1 = (sd_7[1] * (z_7[1]));
  r_7_logvB_1 = r_7;
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_logK1 = X_logK1 * b_logK1;
    // initialize linear predictor term
    vector[N] nlp_logBPnd = X_logBPnd * b_logBPnd;
    // initialize linear predictor term
    vector[N] nlp_logVnd = X_logVnd * b_logVnd;
    // initialize linear predictor term
    vector[N] nlp_logk4 = X_logk4 * b_logk4;
    // initialize linear predictor term
    vector[N] nlp_logvB = X_logvB * b_logvB;
    // initialize non-linear predictor term
    vector[N] mu;
    // initialize linear predictor term
    vector[N] sigma = X_sigma * b_sigma + Xs_sigma * bs_sigma + Zs_sigma_1_1 * s_sigma_1_1;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_logK1[n] += r_3_logK1_1[J_3[n]] * Z_3_logK1_1[n] + r_4_logK1_1[J_4[n]] * Z_4_logK1_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_logBPnd[n] += r_3_logBPnd_2[J_3[n]] * Z_3_logBPnd_2[n] + r_4_logBPnd_2[J_4[n]] * Z_4_logBPnd_2[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_logVnd[n] += r_3_logVnd_3[J_3[n]] * Z_3_logVnd_3[n] + r_4_logVnd_3[J_4[n]] * Z_4_logVnd_3[n] + r_5_logVnd_1[J_5[n]] * Z_5_logVnd_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_logk4[n] += r_3_logk4_4[J_3[n]] * Z_3_logk4_4[n] + r_4_logk4_4[J_4[n]] * Z_4_logk4_4[n] + r_5_logk4_2[J_5[n]] * Z_5_logk4_2[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_logvB[n] += r_6_logvB_1[J_6[n]] * Z_6_logvB_1[n] + r_7_logvB_1[J_7[n]] * Z_7_logvB_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      sigma[n] += r_1_sigma_1[J_1[n]] * Z_1_sigma_1[n] + r_2_sigma_1[J_2[n]] * Z_2_sigma_1[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      sigma[n] = exp(sigma[n]);
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = twotcm_log_stan(nlp_logK1[n] , nlp_logVnd[n] , nlp_logBPnd[n] , nlp_logk4[n] , nlp_logvB[n] , C_1[n] , C_2[n] , C_3[n] , C_4[n] , C_5[n] , C_6[n] , C_7[n] , C_8[n] , C_9[n] , C_10[n] , C_11[n]);
    }
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += normal_lpdf(b_logK1[1] | -2.5, 0.25);
  target += normal_lpdf(b_logK1[2] | 0, 0.3);
  target += normal_lpdf(b_logK1[3] | 0, 0.3);
  target += normal_lpdf(b_logK1[4] | 0, 0.3);
  target += normal_lpdf(b_logK1[5] | 0, 0.3);
  target += normal_lpdf(b_logK1[6] | 0, 0.3);
  target += normal_lpdf(b_logK1[7] | 0, 0.3);
  target += normal_lpdf(b_logK1[8] | 0, 0.3);
  target += normal_lpdf(b_logK1[9] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[1] |  2, 0.25);
  target += normal_lpdf(b_logBPnd[2] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[3] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[4] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[5] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[6] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[7] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[8] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[9] | 0, 0.3);
  target += normal_lpdf(b_logBPnd[10] | 0, 0.2);
  target += normal_lpdf(b_logVnd | -1, 0.25);
  target += normal_lpdf(b_logk4 | -4, 0.25);
  target += normal_lpdf(b_logvB | -4, 0.5);
  target += normal_lpdf(b_sigma | -5, 1);
  target += student_t_lpdf(bs_sigma[1] | 3, 0, 4)
    - 1 * student_t_lccdf(0 | 3, 0, 4);
  target += student_t_lpdf(sds_sigma_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  // target += std_normal_lpdf(zs_sigma_1_1);
  for(n in 1:knots_sigma_1[1]) {
    target += normal_lpdf(rs_sigma_1_1[n] |  0, sds_sigma_1_1);
  }
  target += normal_lpdf(sd_1 | 0, 0.5)
    - 1 * normal_lccdf(0 | 0, 0.5);
  // target += std_normal_lpdf(z_1[1]);
  for(n in 1:N_1) {
    target += normal_lpdf(r_1[n] |  0, sd_1);
  }
  target += normal_lpdf(sd_2 | 0, 0.3)
    - 1 * normal_lccdf(0 | 0, 0.3);
  // target += std_normal_lpdf(z_2[1]);
  for(n in 1:N_2) {
    target += normal_lpdf(r_2[n] |  0, sd_2);
  }
  target += normal_lpdf(sd_3 | 0, 0.3)
    - 4 * normal_lccdf(0 | 0, 0.3);
  // target += std_normal_lpdf(to_vector(z_3));
  for(n in 1:N_3) {
    target += multi_normal_cholesky_lpdf(r_3[n,] |  
                rep_row_vector(0.0, M_3), diag_pre_multiply(sd_3, L_3));
  }
  target += lkj_corr_cholesky_lpdf(L_3 | 1);
  target += normal_lpdf(sd_4 | 0, 0.025)
    - 4 * normal_lccdf(0 | 0, 0.025);
  target += std_normal_lpdf(to_vector(z_4));
  target += lkj_corr_cholesky_lpdf(L_4 | 2);
  target += normal_lpdf(sd_5 | 0, 0.1)
    - 2 * normal_lccdf(0 | 0, 0.1);
  // target += std_normal_lpdf(to_vector(z_5));
  for(n in 1:N_5) {
    target += normal_lpdf(r_5[n] |  0, sd_5);
  }
  target += lkj_corr_cholesky_lpdf(L_5 | 2);
  target += normal_lpdf(sd_6 | 0, 0.1)
    - 1 * normal_lccdf(0 | 0, 0.1);
  target += std_normal_lpdf(z_6[1]);
  // for(n in 1:N_6) {
  //   target += normal_lpdf(r_6[n] |  0, sd_6);
  // }
  target += normal_lpdf(sd_7 | 0, 0.5)
    - 1 * normal_lccdf(0 | 0, 0.5);
  // target += std_normal_lpdf(z_7[1]);
  for(n in 1:N_7) {
    target += normal_lpdf(r_7[n] |  0, sd_7);
  }
}
generated quantities {
  // compute group-level correlations
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  // compute group-level correlations
  corr_matrix[M_4] Cor_4 = multiply_lower_tri_self_transpose(L_4);
  vector<lower=-1,upper=1>[NC_4] cor_4;
  // compute group-level correlations
  corr_matrix[M_5] Cor_5 = multiply_lower_tri_self_transpose(L_5);
  vector<lower=-1,upper=1>[NC_5] cor_5;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_4) {
    for (j in 1:(k - 1)) {
      cor_4[choose(k - 1, 2) + j] = Cor_4[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_5) {
    for (j in 1:(k - 1)) {
      cor_5[choose(k - 1, 2) + j] = Cor_5[j, k];
    }
  }
}

