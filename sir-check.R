gen_sir <- odin2::odin({
  
  # update for the next time step
  
  update(N_1) <- S_1 + I_1 + R_1
  update(S_1) <- S_1 + n_birth_1 - n_SI_event_1
  update(I_1) <- I_1 + n_SI_1 - n_IR_event_1
  update(R_1) <- R_1 + n_IR_1 - n_R_death_1
  
  update(N_2) <- S_2 + I_2 + R_2
  update(S_2) <- S_2 + n_birth_2 - n_SI_event_2
  update(I_2) <- I_2 + n_SI_2 - n_IR_event_2
  update(R_2) <- R_2 + n_IR_2 - n_R_death_2
  
  update(N_3) <- S_3 + I_3 + R_3
  update(S_3) <- S_3 + n_birth_3 - n_SI_event_3
  update(I_3) <- I_3 + n_SI_3 - n_IR_event_3
  update(R_3) <- R_3 + n_IR_3 - n_R_death_3
  
  
  # transition probabilities of (death and I) from S
  
  beta_1 <- beta_0_1 + a1 * cos(2 * pi * time / period) + b1 * sin(2 * pi * time / period) +
    a2 * cos(4 * pi * time / period) + b2 * sin(4 * pi * time / period)
  
  beta_2 <- beta_0_2 + a1 * cos(2 * pi * time / period) + b1 * sin(2 * pi * time / period) +
    a2 * cos(4 * pi * time / period) + b2 * sin(4 * pi * time / period)
  
  beta_3 <- beta_0_3 + a1 * cos(2 * pi * time / period) + b1 * sin(2 * pi * time / period) +
    a2 * cos(4 * pi * time / period) + b2 * sin(4 * pi * time / period)
  
  foi_1 <- if ((beta_1 * (I_1 + w_1) / N_1) < 0) 0 else (beta_1 * (I_1 + w_1) / N_1)
  foi_2 <- if ((beta_2 * (I_2 + w_2) / N_2) < 0) 0 else (beta_2 * (I_2 + w_2) / N_2)
  foi_3 <- if ((beta_3 * (I_3 + w_3) / N_3) < 0) 0 else (beta_3 * (I_3 + w_3) / N_3)
  
  # from S
  p_SI_event_1 <- 1 - exp(-(foi_1 + death_rate) * dt)  
  p_SI_death_1 <- death_rate / (foi_1 + death_rate)
  n_SI_event_1 <- if (S_1 < 1) 0 else Binomial(S_1, p_SI_event_1)
  n_SI_death_1 <- if (n_SI_event_1 < 1) 0 else Binomial(n_SI_event_1, p_SI_death_1)
  n_SI_1 <- n_SI_event_1 - n_SI_death_1
  
  p_SI_event_2 <- 1 - exp(-(foi_2 + death_rate) * dt)  
  p_SI_death_2 <- death_rate / (foi_2 + death_rate)
  n_SI_event_2 <- if (S_2 < 1) 0 else Binomial(S_2, p_SI_event_2)
  n_SI_death_2 <- if (n_SI_event_2 < 1) 0 else Binomial(n_SI_event_2, p_SI_death_2)
  n_SI_2 <- n_SI_event_2 - n_SI_death_2
  
  p_SI_event_3 <- 1 - exp(-(foi_3 + death_rate) * dt)  
  p_SI_death_3 <- death_rate / (foi_3 + death_rate)
  n_SI_event_3 <- if (S_3 < 1) 0 else Binomial(S_3, p_SI_event_3)
  n_SI_death_3 <- if (n_SI_event_3 < 1) 0 else Binomial(n_SI_event_3, p_SI_death_3)
  n_SI_3 <- n_SI_event_3 - n_SI_death_3
  
  # from I
  p_IR_event <- 1 - exp(-(gamma + death_rate) * dt)  
  p_IR_death <- death_rate / (death_rate + gamma)
  
  n_IR_event_1 <- if (I_1 < 1) 0 else Binomial(I_1, p_IR_event)
  n_IR_death_1 <- if (n_IR_event_1 < 1) 0 else Binomial(n_IR_event_1, p_IR_death)
  n_IR_1 <- n_IR_event_1 - n_IR_death_1
  
  n_IR_event_2 <- if (I_2 < 1) 0 else Binomial(I_2, p_IR_event)
  n_IR_death_2 <- if (n_IR_event_2 < 1) 0 else Binomial(n_IR_event_2, p_IR_death)
  n_IR_2 <- n_IR_event_2 - n_IR_death_2
  
  n_IR_event_3 <- if (I_3 < 1) 0 else Binomial(I_3, p_IR_event)
  n_IR_death_3 <- if (n_IR_event_3 < 1) 0 else Binomial(n_IR_event_3, p_IR_death)
  n_IR_3 <- n_IR_event_3 - n_IR_death_3
  
  # from R 
  n_R_death_1 <- if (R_1 < 1) 0 else Binomial(R_1, 1 - exp(-death_rate * dt)) 
  n_R_death_2 <- if (R_2 < 1) 0 else Binomial(R_2, 1 - exp(-death_rate * dt)) 
  n_R_death_3 <- if (R_3 < 1) 0 else Binomial(R_3, 1 - exp(-death_rate * dt)) 
  
  # from N
  n_birth_1 <- if (N_1 < 1) 0 else Binomial(N_1, 1 - exp(-birth_rate * dt))
  n_birth_2 <- if (N_2 < 1) 0 else Binomial(N_2, 1 - exp(-birth_rate * dt))
  n_birth_3 <- if (N_3 < 1) 0 else Binomial(N_3, 1 - exp(-birth_rate * dt))
  
  # Initial values
  
  initial(N_1) <- N_0
  initial(N_2) <- N_0
  initial(N_3) <- N_0
  initial(S_1) <- floor(N_0 * s_0_1)
  initial(S_2) <- floor(N_0 * s_0_2)
  initial(S_3) <- floor(N_0 * s_0_3)
  initial(I_1) <- floor(N_0 * i_0_1)
  initial(I_2) <- floor(N_0 * i_0_2)
  initial(I_3) <- floor(N_0 * i_0_3)
  initial(R_1) <- N_0 - floor(N_0 * s_0_1) - floor(N_0 * i_0_1)
  initial(R_2) <- N_0 - floor(N_0 * s_0_2) - floor(N_0 * i_0_2)
  initial(R_3) <- N_0 - floor(N_0 * s_0_3) - floor(N_0 * i_0_3)
  

  # Parameters unit daily
  
  period <- parameter()     
  gamma <- parameter()      
  N_0 <- parameter()        
  
  beta_0_1 <- parameter()     
  a1 <- parameter()
  a2 <- parameter()
  b1 <- parameter()
  b2 <- parameter()
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
  
  
  p_1 <- parameter()
  p_2 <- parameter()
  p_3 <- parameter()

  noise <- Exponential(rate = 1e6)
  
  update(incidence_1) <- n_SI_1 + incidence_1
  initial(incidence_1, zero_every = 1) <- 0
  
  update(incidence_2) <- n_SI_2 + incidence_2
  initial(incidence_2, zero_every = 1) <- 0
  
  update(incidence_3) <- n_SI_3 + incidence_3
  initial(incidence_3, zero_every = 1) <- 0
  
  update(incidence_tol) <- incidence_1 * p_1 + incidence_2 * p_2 + incidence_3 * p_3 + noise
  initial(incidence_tol) <- 0
  
  # 1. error !!!
  size <- parameter()
  update(inc_cases) <- NegativeBinomial(size = size, mu = incidence_tol)
  initial(inc_cases) <- 0
  
  # update(inc_cases) <- Poisson(incidence_tol)
  # initial(inc_cases) <- 0
  
  # 2. multi-stream data likelihood !!!
  
  # observation model

  cases_hfmd <- data()
  cases_hfmd ~ NegativeBinomial(size = size, mu = incidence_tol)
  # cases_hfmd ~ Poisson(incidence_tol)
  
  rho_1 <- incidence_1*p_1 / incidence_tol
  rho_2 <- incidence_2*p_2 / incidence_tol
  rho_3 <- incidence_3*p_3 / incidence_tol

  cases_pos <- data()

  cases_1 <- data()
  cases_2 <- data()
  cases_3 <- data()

  cases_1 ~ Binomial(cases_pos, rho_1)
  cases_2 ~ Binomial(cases_pos, rho_2)
  cases_3 ~ Binomial(cases_pos, rho_3)
  
})

