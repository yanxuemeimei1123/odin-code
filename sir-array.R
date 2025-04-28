gen_sir <- odin2::odin({
  
  # update for the next time step
  
  update(N[]) <- S[i] + I[i] + R[i]
  update(S[]) <- S[i] + n_birth[i] - n_SI_event[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR_event[i]
  update(R[]) <- R[i] + n_IR[i] - n_R_death[i]
  
  
  # transition probabilities of (death and I) from S
  
  beta[1] <- beta_0_1 + a1 * cos(2 * pi * time / period) + b1 * sin(2 * pi * time / period) +
    a2 * cos(4 * pi * time / period) + b2 * sin(4 * pi * time / period)
  
  beta[2] <- beta_0_2 + a1 * cos(2 * pi * time / period) + b1 * sin(2 * pi * time / period) +
    a2 * cos(4 * pi * time / period) + b2 * sin(4 * pi * time / period)
  
  beta[3] <- beta_0_3 + a1 * cos(2 * pi * time / period) + b1 * sin(2 * pi * time / period) +
    a2 * cos(4 * pi * time / period) + b2 * sin(4 * pi * time / period)
  
  foi[1] <- beta[1] * (I[1] + w_1) / N[1]
  foi[2] <- beta[2] * (I[2] + w_2) / N[2]
  foi[3] <- beta[3] * (I[3] + w_3) / N[3]
  
  # from S
  p_SI_event[] <- 1 - exp(-(foi[i] + death_rate) * dt)  
  p_SI_death[] <- death_rate / (foi[i] + death_rate)
  
  n_SI_event[] <- Binomial(S[i], p_SI_event[i])
  n_SI_death[] <- Binomial(n_SI_event[i], p_SI_death[i])
  n_SI[] <- n_SI_event[i] - n_SI_death[i]
  
  
  # from I
  p_IR_event <- 1 - exp(-(gamma + death_rate) * dt)  
  p_IR_death <- death_rate / (death_rate + gamma)
  
  n_IR_event[] <- Binomial(I[i], p_IR_event)
  n_IR_death[] <- Binomial(n_IR_event[i], p_IR_death)
  n_IR[] <- n_IR_event[i] - n_IR_death[i]
  
  
  # from R 
  n_R_death[] <- Binomial(R[i], 1 - exp(-death_rate * dt)) 
  
  
  # from N
  n_birth[] <- Binomial(N[i], 1 - exp(-birth_rate * dt))
  
  
  # Initial values
  
  initial(N[]) <- N_0
  initial(S[1]) <- floor(N_0 * s_0_1)
  initial(S[2]) <- floor(N_0 * s_0_2)
  initial(S[3]) <- floor(N_0 * s_0_3)
  
  initial(I[1]) <- floor(N_0 * i_0_1)
  initial(I[2]) <- floor(N_0 * i_0_2)
  initial(I[3]) <- floor(N_0 * i_0_3)
  
  initial(R[1]) <- N_0 - floor(N_0 * s_0_1) - floor(N_0 * i_0_1)
  initial(R[2]) <- N_0 - floor(N_0 * s_0_2) - floor(N_0 * i_0_2)
  initial(R[3]) <- N_0 - floor(N_0 * s_0_3) - floor(N_0 * i_0_3)
  
  
  # Parameters unit daily
  
  period <- parameter()     
  gamma <- parameter()      
  N_0 <- parameter()        
  a1 <- parameter()
  a2 <- parameter()
  b1 <- parameter()
  b2 <- parameter()
  
  beta_0_1 <- parameter()     
  s_0_1 <- parameter()        
  i_0_1 <- parameter()        
  w_1 <- parameter()          
  
  beta_0_2 <- parameter()     
  s_0_2 <- parameter()        
  i_0_2 <- parameter()        
  w_2 <- parameter() 
  
  beta_0_3 <- parameter()     
  s_0_3 <- parameter()        
  i_0_3 <- parameter()        
  w_3 <- parameter() 
  
  birth_rate <- interpolate(birth_time, birth_value, "spline")
  birth_time <- parameter(constant = TRUE)
  birth_value <- parameter(constant = TRUE)
  dim(birth_time, birth_value) <- parameter(rank = 1)
  
  death_rate <- interpolate(death_time, death_value, "spline")
  death_time <- parameter(constant = TRUE)
  death_value <- parameter(constant = TRUE)
  dim(death_time, death_value) <- parameter(rank = 1)
  
  
  # Dimensions of arrays
  
  n_sero <- parameter(3)
  
  dim(N) <- n_sero
  dim(S, n_birth, n_SI_event, n_SI_death, n_SI, p_SI_event, p_SI_death) <- n_sero
  dim(I, n_IR_event, n_IR_death, n_IR) <- n_sero
  dim(R, n_R_death) <- n_sero
  dim(beta, foi) <- n_sero
  dim(incidence) <- n_sero
  
  
  # observation model

  p_1 <- parameter()
  p_2 <- parameter()
  p_3 <- parameter()
  size <- parameter()
  noise <- Exponential(rate = 1e6)
  
  update(incidence[]) <- n_SI[i] + incidence[i]
  initial(incidence[]) <- 0
  # initial(incidence[], zero_every = 1) <- 0  # error !!!
  
  update(incidence_tol) <- incidence[1] * p_1 + incidence[2] * p_2 + incidence[3] * p_3 + noise
  initial(incidence_tol) <- 0
  
  # update(inc_cases) <- NegativeBinomial(size = size, mu = incidence_tol)
  # initial(inc_cases) <- 0

  # observation model

  cases_hfmd <- data()
  cases_hfmd ~ NegativeBinomial(size = size, mu = incidence_tol)

  rho_1 <- incidence[1] * p_1 / incidence_tol
  rho_2 <- incidence[2] * p_2 / incidence_tol
  rho_3 <- incidence[3] * p_3 / incidence_tol

  cases_1 <- data()
  cases_2 <- data()
  cases_3 <- data()
  cases_pos <- data()

  cases_1 ~ Binomial(cases_pos, rho_1)
  cases_2 ~ Binomial(cases_pos, rho_2)
  cases_3 ~ Binomial(cases_pos, rho_3)

})

