
pacman::p_load(tidyverse, 
               deSolve, 
               ggplot2, 
               glue, 
               readxl, 
               data.table, 
               lubridate, 
               Rcpp, 
               rlang
)

cppFunction('
#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector EQ_omp(const NumericVector transit,
                     const IntegerVector iu1,
                     const IntegerVector iu2,
                     const IntegerVector iv1,
                     const IntegerVector iv2,
                     const int neq) {

  const int n = transit.size();
  int nthreads = 1;
  #ifdef _OPENMP
    nthreads = omp_get_max_threads();
  #endif

  // thread-local accumulators: nthreads x neq
  std::vector<std::vector<double>> local(nthreads, std::vector<double>(neq, 0.0));

  #pragma omp parallel
  {
    int tid = 0;
    #ifdef _OPENMP
      tid = omp_get_thread_num();
    #endif
    std::vector<double> &acc = local[tid];

    #pragma omp for schedule(static)
    for (int i = 0; i < n; ++i) {
      const int j1 = iu1[i] - 1;
      const int j2 = iu2[i] - 1;
      const double v = transit[i];
      acc[j1] += v * iv1[i];
      acc[j2] += v * iv2[i];
    }
  }

  // reduce to a single output vector
  NumericVector out(neq);
  for (int t = 0; t < nthreads; ++t)
    for (int k = 0; k < neq; ++k)
      out[k] += local[t][k];

  return out;
}
')


N = 12
B = 14 + 8*3 + 3  # number of variables per patch
A = 30  + # Transitions for Child
  15*3 + # Transitions for Men, Non Pregnant Women and Pregnant Women
  ((14 + 8*3)*(N)) + # Migration 
  (14 + 8*3) + # Deaths
  (8 + 8 + 8*2) + # Aging for children and women moving from not pregnant to pregnant 
  (3*4) # Additional Deaths from complicated, uncomplicated and asymptomatic infection
3 + # Mosquito movements
  3 # Mosquito Deaths
V = N*B
L = N*A

startyear = 2022 # starting year of simulation 2022-01-01

tyears = 15 # total Years of simulation 
dtout = 1/52 # output timestep
tsteps = round(tyears/dtout) # number of time steps
time = startyear+seq(0, tyears, dtout) # time vector

alivepopHum = c(1:(B - 3))
alivepopMos = c((B - 2):B)
popc <- rep(0,N)
popM <- rep(0,N)

childClass = c("ch")
otherClasses = c("m","wnp","wp")
compartmentsChildren = c("S", "It", "V1", "V2", "V3", "v4", "vstar", "E", "C", "U", "A", "Tc", "Tu", "R") 
compartmentsOther = c("S", "E", "C", "U", "A", "Tc", "Tu", "R")
var_namesChildren = as.vector(sapply(childClass, function(cl) paste0(compartmentsChildren, "_", cl)))
var_namesOther = as.vector(sapply(otherClasses, function(cl) paste0(compartmentsOther, "_", cl)))

var_namesMosquitos = c("S_Mo", "E_Mo", "I_Mo")

var_names = c(var_namesChildren, var_namesOther, var_namesMosquitos)

getStates = function(x) {
  for (n in 1:N){
    popc[n] <- sum(x[varind[alivepopHum,n]])  
    popM[n] <- sum(x[varind[alivepopMos,n]])  
  }
  
  states = list()
  
  for (j in 1:length(var_names))
  {
    states[[var_names[j]]] = array(x[varind[j, ]], dim = c(N))
  }
  
  states$.pop = c(popc)
  states$.popMos = c(popM)
  
  return(states)
}


varind=matrix(0,nrow=B,ncol=N)
traind=matrix(0,nrow=A,ncol=N)
for (n in 1:N){
  for (b in 1:B){
    varind[b,n]<-(n-1)*B+b
  }
  for (a in 1:A){
    traind[a,n]<-(n-1)*A+a
  }
}

migrationVector = c(1:N)

transitions = matrix(0, nrow=L, ncol=4)
for (n in 1:N) {
  indAdd = 1
  
  # -> It
  transitions[traind[indAdd,n],] <- c(varind[1,n],  0, varind[2,n], +1) 
  indAdd = indAdd +1
  
  # It -> S
  transitions[traind[indAdd,n],] <- c(varind[2,n],  -1, varind[1,n], +1) 
  indAdd = indAdd + 1
  
  # It -> V1 
  transitions[traind[indAdd,n],] <- c(varind[2,n],  -1, varind[3,n], +1) 
  indAdd = indAdd + 1
  # V1 -> V2 
  transitions[traind[indAdd,n],] <- c(varind[3,n],  -1, varind[4,n], +1) 
  indAdd = indAdd + 1
  # V2 -> V3 
  transitions[traind[indAdd,n],] <- c(varind[4,n],  -1, varind[5,n], +1) 
  indAdd = indAdd + 1
  
  # V3 -> V4 
  transitions[traind[indAdd,n],] <- c(varind[5,n],  -1, varind[6,n], +1) 
  indAdd = indAdd + 1
  
  # V4 -> V* 
  transitions[traind[indAdd,n],] <- c(varind[6,n],  -1, varind[7,n], +1) 
  indAdd = indAdd + 1
  
  
  # V1 -> E
  transitions[traind[indAdd,n],] <- c(varind[3,n],  -1, varind[8,n], +1) 
  indAdd = indAdd + 1
  # V2 -> E
  transitions[traind[indAdd,n],] <- c(varind[4,n],  -1, varind[8,n], +1) 
  indAdd = indAdd + 1
  # V3 -> E
  transitions[traind[indAdd,n],] <- c(varind[5,n],  -1, varind[8,n], +1) 
  indAdd = indAdd + 1
  # V4 -> E
  transitions[traind[indAdd,n],] <- c(varind[6,n],  -1, varind[8,n], +1) 
  indAdd = indAdd + 1
  # V* -> E
  transitions[traind[indAdd,n],] <- c(varind[7,n],  -1, varind[8,n], +1) 
  indAdd = indAdd + 1
  
  # V1 -> S
  transitions[traind[indAdd,n],] <- c(varind[3,n],  -1, varind[1,n], +1)
  indAdd = indAdd + 1
  # V2 -> S
  transitions[traind[indAdd,n],] <- c(varind[4,n],  -1, varind[1,n], +1) 
  indAdd = indAdd + 1
  # V3 -> V*
  transitions[traind[indAdd,n],] <- c(varind[5,n],  -1, varind[8,n], +1) 
  indAdd = indAdd + 1
  # V* -> S
  transitions[traind[indAdd,n],] <- c(varind[7,n],  -1, varind[1,n], +1) 
  indAdd = indAdd + 1
  
  # Death
  for (i in 1:B)
  {
    transitions[traind[indAdd + i,n],] <- c(varind[i,n],  -1, varind[1,n], 0) #  X -> 
  }
  indAdd = indAdd + B
  print(indAdd)
  
  
  
  # Infection pathway  
  
  for (j in 1:4)
  {
    
    if (j == 1)
    {
      
      transitions[traind[1 + indAdd, n],] <- c(varind[1,n], -1, varind[8,n], +1) # S -> E
      transitions[traind[2 + indAdd, n],] <- c(varind[8,n], -1, varind[9,n], +1) # E -> C
      transitions[traind[3 + indAdd, n],] <- c(varind[8,n], -1, varind[10,n], +1) # E -> U
      transitions[traind[4 + indAdd, n],] <- c(varind[8,n], -1, varind[11,n], +1) # E -> Tc
      transitions[traind[3 + indAdd, n],] <- c(varind[8,n], -1, varind[12,n], +1) # E -> Tu
      transitions[traind[4 + indAdd, n],] <- c(varind[8,n], -1, varind[13,n], +1) # E -> A
      
      transitions[traind[4 + indAdd, n],] <- c(varind[9,n], -1, varind[13,n], +1) # C -> A
      
      transitions[traind[4 + indAdd, n],] <- c(varind[10,n], -1, varind[13,n], +1) # U -> A
      
      transitions[traind[4 + indAdd, n],] <- c(varind[11,n], -1, varind[13,n], +1) # Tc -> A
      transitions[traind[4 + indAdd, n],] <- c(varind[11,n], -1, varind[14,n], +1) # Tc -> R
      
      transitions[traind[4 + indAdd, n],] <- c(varind[12,n], -1, varind[13,n], +1) # Tu -> A
      transitions[traind[4 + indAdd, n],] <- c(varind[12,n], -1, varind[14,n], +1) # Tu -> R
      
      transitions[traind[4 + indAdd, n],] <- c(varind[13,n], -1, varind[14,n], +1) # A -> R
      
      transitions[traind[4 + indAdd, n],] <- c(varind[14,n], -1, varind[1,n], +1) # R -> S
      
    
      
      indAdd = indAdd + 14
      
    }
    
    if (j > 1)
    {
      transitions[traind[1 + indAdd, n],] <- c(varind[1,n], -1, varind[2,n], +1) # S -> E
      transitions[traind[2 + indAdd, n],] <- c(varind[2,n], -1, varind[3,n], +1) # E -> C
      transitions[traind[3 + indAdd, n],] <- c(varind[2,n], -1, varind[4,n], +1) # E -> U
      transitions[traind[4 + indAdd, n],] <- c(varind[2,n], -1, varind[5,n], +1) # E -> Tc
      transitions[traind[3 + indAdd, n],] <- c(varind[2,n], -1, varind[6,n], +1) # E -> Tu
      transitions[traind[4 + indAdd, n],] <- c(varind[2,n], -1, varind[7,n], +1) # E -> A
      
      transitions[traind[4 + indAdd, n],] <- c(varind[3,n], -1, varind[7,n], +1) # C -> A
      
      transitions[traind[4 + indAdd, n],] <- c(varind[4,n], -1, varind[7,n], +1) # U -> A
      
      transitions[traind[4 + indAdd, n],] <- c(varind[5,n], -1, varind[7,n], +1) # Tc -> A
      transitions[traind[4 + indAdd, n],] <- c(varind[5,n], -1, varind[8,n], +1) # Tc -> R
      
      transitions[traind[4 + indAdd, n],] <- c(varind[6,n], -1, varind[7,n], +1) # Tu -> A
      transitions[traind[4 + indAdd, n],] <- c(varind[6,n], -1, varind[8,n], +1) # Tu -> R
      
      transitions[traind[4 + indAdd, n],] <- c(varind[7,n], -1, varind[8,n], +1) # A -> R
      
      transitions[traind[4 + indAdd, n],] <- c(varind[8,n], -1, varind[1,n], +1) # R -> S
      
      indAdd = indAdd + 14
      
    }
    
  }
  
  # Aging child to man
  
  transitions[traind[1 + indAdd, n],] <- c(varind[1,n], -1, varind[15,n], +1)
  transitions[traind[2 + indAdd, n],] <- c(varind[8,n], -1, varind[16,n], +1)
  transitions[traind[3 + indAdd, n],] <- c(varind[9,n], -1, varind[17,n], +1)
  transitions[traind[4 + indAdd, n],] <- c(varind[10,n], -1, varind[18,n], +1)
  transitions[traind[5 + indAdd, n],] <- c(varind[11,n], -1, varind[19,n], +1)
  transitions[traind[6 + indAdd, n],] <- c(varind[12,n], -1, varind[20,n], +1)
  transitions[traind[7 + indAdd, n],] <- c(varind[13,n], -1, varind[21,n], +1)
  transitions[traind[8 + indAdd, n],] <- c(varind[14,n], -1, varind[22,n], +1)

  
  indAdd = indAdd + 8
  
  # Aging child to woman
  
  transitions[traind[1 + indAdd, n],] <- c(varind[1,n], -1, varind[27,n], +1)
  transitions[traind[2 + indAdd, n],] <- c(varind[8,n], -1, varind[28,n], +1)
  transitions[traind[3 + indAdd, n],] <- c(varind[9,n], -1, varind[29,n], +1)
  transitions[traind[4 + indAdd, n],] <- c(varind[10,n], -1, varind[30,n], +1)
  transitions[traind[5 + indAdd, n],] <- c(varind[11,n], -1, varind[31,n], +1)
  transitions[traind[6 + indAdd, n],] <- c(varind[12,n], -1, varind[32,n], +1)
  transitions[traind[7 + indAdd, n],] <- c(varind[13,n], -1, varind[33,n], +1)
  transitions[traind[8 + indAdd, n],] <- c(varind[14,n], -1, varind[34,n], +1)

  
  indAdd = indAdd + 8
  
  
  
  # Pregnancy
  
  for (i in 1:8)
  {
    transitions[traind[i + indAdd,n],] <- c(varind[22 + i,n],  -1, varind[30 + i,n], +1) # Not Pregnant to Pregnant
    
  }
  indAdd = indAdd + 8
  
  
  for (i in 1:8)
  {
    transitions[traind[i + indAdd,n],] <- c(varind[30 + i,n],  -1, varind[22 + i,n], +1) # Pregnant to Not Pregnant
  }
  
  indAdd = indAdd + 8
  
  
  
  # Migration
  
  k = 1
  for (z in 1:N)
  {
    for (i in 1:(B-3))
    {
      transitions[traind[indAdd + 1,n],] <- c(varind[i,n],  -1, varind[i,migrationVector[z]], +1) #  Migrating from n to every 1 to 12, excluding n
      #print(indAdd + k)
      indAdd = indAdd + 1
    }
    
  }
  
  # Additional Mortality
  
  # Children
  transitions[traind[1 + indAdd, n],] <- c(varind[9,n], -1, varind[1,n], 0)
  transitions[traind[2 + indAdd, n],] <- c(varind[10,n], -1, varind[1,n], 0)
  transitions[traind[3 + indAdd, n],] <- c(varind[11,n], -1, varind[1,n], 0)
  transitions[traind[4 + indAdd, n],] <- c(varind[12,n], -1, varind[1,n], 0)
  transitions[traind[5 + indAdd, n],] <- c(varind[13,n], -1, varind[1,n], 0)
  
  # Men
  transitions[traind[6 + indAdd, n],] <- c(varind[17,n], -1, varind[1,n], 0)
  transitions[traind[7 + indAdd, n],] <- c(varind[18,n], -1, varind[1,n], 0)
  transitions[traind[8 + indAdd, n],] <- c(varind[19,n], -1, varind[1,n], 0)
  transitions[traind[9 + indAdd, n],] <- c(varind[20,n], -1, varind[1,n], 0)
  transitions[traind[10 + indAdd, n],] <- c(varind[21,n], -1, varind[1,n], 0)
  
  
  # Women Not Pregnant
  transitions[traind[11 + indAdd, n],] <- c(varind[25,n], -1, varind[1,n], 0)
  transitions[traind[12 + indAdd, n],] <- c(varind[26,n], -1, varind[1,n], 0)
  transitions[traind[13 + indAdd, n],] <- c(varind[27,n], -1, varind[1,n], 0)
  transitions[traind[14 + indAdd, n],] <- c(varind[28,n], -1, varind[1,n], 0)
  transitions[traind[15 + indAdd, n],] <- c(varind[29,n], -1, varind[1,n], 0)
  
  # Women Pregnant
  transitions[traind[16 + indAdd, n],] <- c(varind[33,n], -1, varind[1,n], 0)
  transitions[traind[17 + indAdd, n],] <- c(varind[34,n], -1, varind[1,n], 0)
  transitions[traind[18 + indAdd, n],] <- c(varind[35,n], -1, varind[1,n], 0)
  transitions[traind[19 + indAdd, n],] <- c(varind[36,n], -1, varind[1,n], 0)
  transitions[traind[20 + indAdd, n],] <- c(varind[37,n], -1, varind[1,n], 0)
  
  # Mosquito
  
  transitions[traind[21 + indAdd,n],] <- c(varind[1,n],  0, varind[39,n], +1) # Birth
  # Mosquito deaths already included above
  
  transitions[traind[22 + indAdd,n],] <- c(varind[39,n],  -1, varind[40,n], 1) # S -> E
  transitions[traind[23 + indAdd,n],] <- c(varind[40,n],  -1, varind[41,n], 1) # E -> I
  
  
  
}

transitionsiu1 <- transitions[,1]
transitionsiv1 <- transitions[,2]
transitionsiu2 <- transitions[,3]
transitionsiv2 <- transitions[,4]

inputs<-function(parameters, scenario){
  
}

parsHuman = readxl::read_excel("SBRMIC008_MMID_App/Parameters.xlsx", sheet = "Humans") %>% 
  pull(var=Value, name=Name)
migrationMatrix = readxl::read_excel("SBRMIC008_MMID_App/Parameters.xlsx", sheet = "migrationMatrix")
migrationMatrix = migrationMatrix[, 2:ncol(migrationMatrix)]
parsMosquito = readxl::read_excel("SBRMIC008_MMID_App/Parameters.xlsx", sheet = "Mosquitos") %>% 
  pull(var=Value, name=Name)

parameters = c(parsHuman, parsMosquito)



get_comp_vecs <- function(names, suffix, env = parent.frame()) {
  sapply(names, function(nm) rlang::eval_tidy(rlang::sym(paste0(nm, "_", suffix)), env = env))
}

malrates <- function(x, input, parameters, distMat, t, ti, scenario) {
  
  pars <- as.list(parameters)
  
  
  pars[names(scenario)] <- scenario
  
  dat <- c(pars, getStates(x))
  #data_for_with <- c(as.list(parameters), getStates(x), scenario)
  with(dat, {
    
    totalAsymp = A_ch + A_m + A_wnp + A_wp
    totalTreatU = Tu_ch + Tu_m + Tu_wnp + Tu_wp
    totalTreatC = Tc_ch + Tc_m + Tc_wnp + Tc_wp
    totalU = U_ch + U_m + U_wnp + U_wp
    totalC = C_ch + C_m + C_wnp + C_wp
    
    totalS = S_ch + S_m + S_wnp + S_wp
    totalE = E_ch + E_m + E_wnp + E_wp
    totalR = R_ch + R_m + R_wnp + R_wp
    
    totalImmuneAndVaccine = It_ch + v1_ch + v2_ch + v3_ch + v4_ch + vstar_ch
    
    N_Hum = totalS + totalE + totalR + totalImmuneAndVaccine + totalAsymp + totalTreatU + totalTreatC + totalU + totalC
    
    totalInfected = totalAsymp + totalTreatU + totalTreatC + totalU + totalC
    
    lambda_hum = beta_Mos*(I_Mo/N_Mo)*beta_BN*pBN
    
    lambda_mos = beta_Mos*((totalInfected)/N_Hum)*Beta_BN*pBN
    
    b_mos = b_mos*beta_Temp*beta_Rain*beta_BN_b*beta_Spray 
    
    N_Mos = S_Mo + E_Mo + I_Mo
    
    print("#############################")
    print(lambda_v)
    print(lambda_mos_V)
    
    
    #notPregnant = 1/40
    # m1 = 0.5
    # m2 = 0.25
    # m3 = 1
    # m4 = 0
    # m5 = 0
    # m6 = 0
    # m7 = 0
    # m8 = 0
    # m9 = 0
    # m10 = 0
    # m11 = 0
    # 
    # m1_m = 0
    # m2_m = 0
    # m3_m = 0
    # m4_m = 0
    # m5_m = 0
    # m6_m = 0
    # m7_m = 0
    # m8_m = 0
    # m9_m = 0
    # m10_m = 0
    # m11_m = 0
    # 
    # migrationFor_ch_wp_wnp = c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)
    # migrationFor_m = c(m1_m, m2_m, m3_m, m4_m, m5_m, m6_m, m7_m, m8_m, m9_m, m10_m, m11_m)
    migrationMatrix = distMat
    
    
    emptyCat_ch = c(S_ch,  It_ch,  v1_ch,  v2_ch,  v3_ch,  v4_ch,  vstar_ch,  E_ch,  C_ch,  U_ch,  A_ch,  Tc_ch,  Tu_ch,  R_ch)
    emptyCat_m = c(S_m,   E_m,  C_m,  U_m,  A_m,  Tc_m,  Tu_m,  R_m)
    emptyCat_wp = c(S_wp,   E_wp,  C_wp,  U_wp,  A_wp,  Tc_wp,  Tu_wp,  R_wp)
    emptyCat_wnp = c(S_wnp,   E_wnp,  C_wnp,  U_wnp,  A_wnp,  Tc_wnp,  Tu_wnp,  R_wnp)
    
    #emptyCat = c(S, Ev, Av, Iv_RDT_TP, Iv_RDT_FN,
    #                Iv_M_TP, Iv_M_FN, Tv, notTv, Tv_D, Ev_D, Rv, 
    #                 EF, AF, IF_RDT_TP, IF_RDT_FN, IF_M_TP, 
    #               IF_M_FN, TF, notTF, RF)
    
    # compartmentsChildren = c("IT", "S", "Ev", "Av", "Iv_RDT_TP", "Iv_RDT_FN", "Iv_M_TP", "Iv_M_FN", 
    #   "Tv", "notTv", "Tv_D", "Ev_D", "Rv", "EF", "AF", "IF_RDT_TP", "IF_RDT_FN", "IF_M_TP", 
    #   "IF_M_FN", "TF", "notTF", "RF")    
    
    # var_namesMosquitos = c("S_Mo", "Ev_Mo", "Iv_Mo", "EF_Mo", "IF_Mo")
    birth = c(b*N_Hum)
    loss_immune = c(eta*IT_ch*(1-pV1))
    
    vaccination = c(eta*It_ch*pV1, delta1*pDelta1*V1_ch, delta2*pDelta2*v2_ch, delta3*pDelta3*v3_ch, delta4*v4_ch)
    vaccineToExposed = c(lambda_hum*beta_child*v1_ch, lambda_hum*beta_child*v2_ch, lambda_hum*beta_child*betaV3*v1_ch, lambda_hum*beta_child*betav4*v1_ch, lambda_hum*beta_child*betaVstar*v1_ch)
    vaccineToSusceptible = c(delta1*(1-pDelta1)*V1_ch, delta2*(1-pDelta2)*v2_ch, delta3*(1-pDelta3)*v3_ch, delta5*vstar_ch)
    
    
    death_ch = mu*emptyCat_ch
    death_m = mu*emptyCat_m
    death_wp = mu*emptyCat_wp
    death_wnp = mu*emptyCat_wnp
    
    
    deathMos = c(mu_mos*S_Mo, mu_mos*Ev_Mo, mu_mos*Iv_Mo, mu_mos*EF_Mo, mu_mos*IF_Mo)
    
    Path_ch = c(lambda_hum*beta_Child*S_ch, sigma*pC*E_ch, sigma*pU*E_ch, sigma*pTC*E_ch, sigma*pTU*E_ch, sigma*pA*E_ch, tauC*C_ch, tauU*U_ch, rC*pNotTC*Tc_ch, rC*(1-pNotTC)*Tc_ch, 
                     rU*pNotTU*Tu_ch, rU*(1-pNotTU)*Tu_ch, delta*A, r*R)
    
    Path_m = c(lambda_hum*S_ch, sigma*pC*E_ch, sigma*pU*E_ch, sigma*pTC*E_ch, sigma*pTU*E_ch, sigma*pA*E_ch, tauC*C_ch, tauU*U_ch, rC*pNotTC*Tc_ch, rC*(1-pNotTC)*Tc_ch, 
                    rU*pNotTU*Tu_ch, rU*(1-pNotTU)*Tu_ch, delta*A, r*R)
    
    Path_wp = c(lambda_hum*beta_Pregnant*S_ch, sigma*pC*E_ch, sigma*pU*E_ch, sigma*pTC*E_ch, sigma*pTU*E_ch, sigma*pA*E_ch, tauC*C_ch, tauU*U_ch, rC*pNotTC*Tc_ch, rC*(1-pNotTC)*Tc_ch, 
                     rU*pNotTU*Tu_ch, rU*(1-pNotTU)*Tu_ch, delta*A, r*R)
    
    Path_wnp = c(lambda_hum*S_ch, sigma*pC*E_ch, sigma*pU*E_ch, sigma*pTC*E_ch, sigma*pTU*E_ch, sigma*pA*E_ch, tauC*C_ch, tauU*U_ch, rC*pNotTC*Tc_ch, rC*(1-pNotTC)*Tc_ch, 
                      rU*pNotTU*Tu_ch, rU*(1-pNotTU)*Tu_ch, delta*A, r*R)
    

    agingToMen = c(a*pMal*S_ch, a*pMal*E_ch, a*pMal*C_ch, a*pMal*U_ch, a*pMal*A_ch, a*pMal*Tc_ch, a*pMal*Tu_ch, a*pMal*R_ch)
    agingToWomen = c(a*pMal*S_ch, a*pMal*E_ch, a*pMal*C_ch, a*pMal*U_ch, a*pMal*A_ch, a*pMal*Tc_ch, a*pMal*Tu_ch, a*pMal*R_ch)

    
    wnpTowp = pregnant*emptyCat_wnp
    wpTownp = notPregnant*emptyCat_wp
    
    
    migrationVec = numeric(N*N*85)
    
    replacementIndex = 1
    
    M_ch  <- matrix(emptyCat_ch,  nrow = N, byrow = FALSE)  # N x 22
    M_m   <- matrix(emptyCat_m,   nrow = N, byrow = FALSE)  # N x 21
    M_wp  <- matrix(emptyCat_wp,  nrow = N, byrow = FALSE)  # N x 21
    M_wnp <- matrix(emptyCat_wnp, nrow = N, byrow = FALSE)  # N x 21
    M_hum <- cbind(M_ch, M_m, M_wp, M_wnp)                  # N x 85  (rows = source patches)
    
    # For each destination z, scale each ROW i by migrationMatrix[i, z]
    blocks <- lapply(1:N, function(z) sweep(M_hum, 1, migrationMatrix[, z], `*`))  # each N x 85
    migration_block <- do.call(cbind, blocks)  # N x (85*N), column blocks in order z = 1..N
    
    # Flatten in the same layout used later (array(..., dim=c(N, A)); c(t(.)))
    migrationVec <- c(t(migration_block))
    
    
    additionalMortCompViv_ch = c(mu_C*C_ch, mu_U*U_ch, mu_Tc*Tc_ch, mu_Tu*Tu_ch, mu_A*A_ch)

    
    additionalMortCompViv_m = c(mu_C*C_m, mu_U*U_m, mu_Tc*Tc_m, mu_Tu*Tu_m, mu_A*A_m)

    
    additionalMortCompViv_wp = c(mu_C*C_wp, mu_U*U_wp, mu_Tc*Tc_wp, mu_Tu*Tu_wp, mu_A*A_wp)

    
    additionalMortCompViv_wnp = c(mu_C*C_wnp, mu_U*U_wnp, mu_Tc*Tc_wnp, mu_Tu*Tu_wnp, mu_A*A_wnp)

    
    # var_namesMosquitos = c("S_Mo", "Ev_Mo", "Iv_Mo", "EF_Mo", "IF_Mo")
    # Mosquito
    
    mosMovement = c(b_mos*N_Mos, lambda_mos*S_Mo, gamma_mos*E_Mo)
    
    tranrate <- array(c(
      birth, loss_immune, vaccination, vaccineToExposed, vaccineToSusceptible, 
      death_ch, death_m, death_wp, death_wnp, deathMos, Path_ch, Path_m, Path_wp, Path_wnp, 
      agingToMen, agingToWomen, wnpTowp, wpTownp, migrationVec, additionalMortCompViv_ch,
      additionalMortCompViv_m, additionalMortCompViv_wp, additionalMortCompViv_wnp,   
      mosMovement), dim=c(N, A))
    tranrate <- c(t(tranrate))
    return(tranrate)
  })
}

postproc <- function(parameters,out,tran) {
  with(as.list(c(parameters)),
       {
         # ************************************************************************************* #
         # for outputting the  time series for each patch
         # ************************************************************************************* #
         
         # Case outputs
         inc_pred<-matrix(0,nrow=length(out[,1]),ncol=N)
         for (n in 1:N){
           inc_pred[,n]<-tran[,c(traind[6,n])]/365  # Incidence
           
         }
         
         
         return(cbind(inc_pred #1
                      
         ))
         
       })
}

# OLD dZ
# dZ <- rep(0.0, V) %>% as.numeric()
# epiModel <- function(t, state, parameters, input, scenario) {
#   with(as.list(c(state, parameters)), {
#     # ************************************************************************************* #
#     # define variables
#     # ************************************************************************************* #
#     
#     Z <- state
#     
#     # rates of change
#     ti <- 1
#     transit <- malrates(state, input, parameters, as.matrix(migrationMatrix), t, ti, scenario)
#     if (any(!is.finite(transit))) {
#       bad <- which(!is.finite(transit))[1]
#       stop(sprintf("non-finite transit at t=%.6f (index %d)", t, bad), call. = FALSE)
#     }
#     
#     EQ(dZ, transit, transitionsiu1, transitionsiu2, transitionsiv1, transitionsiv2)
#     
#     if (any(!is.finite(dZ))) {
#       bad <- which(!is.finite(dZ))[1]
#       stop(sprintf("non-finite dZ at t=%.6f (index %d)", t, bad), call. = FALSE)
#     }
#     
#     list(c(dZ))
#     # return the rate of change
#     # print(c(t, dZ))
#     # browser()
#     list(c(dZ))
#   }
#   ) 
#   # end with(as.list ...
# }

dZ <- rep(0.0, V) %>% as.numeric()
epiModel <- function(t, state, parameters, input, scenario) {
  with(as.list(c(state, parameters)), {
    # ************************************************************************************* #
    # define variables
    # ************************************************************************************* #
    
    Z <- state
    
    # rates of change
    ti <- 1
    transit <- malrates(state, input, parameters, as.matrix(migrationMatrix), t, ti, scenario)
    if (any(!is.finite(transit))) {
      bad <- which(!is.finite(transit))[1]
      stop(sprintf("non-finite transit at t=%.6f (index %d)", t, bad), call. = FALSE)
    }
    
    neq <- length(state)  # or use V
    dZ  <- EQ_omp(transit,
                  transitionsiu1,
                  transitionsiu2,
                  transitionsiv1,
                  transitionsiv2,
                  neq)
    
    if (any(!is.finite(dZ))) {
      bad <- which(!is.finite(dZ))[1]
      stop(sprintf("non-finite dZ at t=%.6f (index %d)", t, bad), call. = FALSE)
    }
    list(dZ)
    # return the rate of change
    # print(c(t, dZ))
    # browser()
    list(c(dZ))
  }
  ) 
  # end with(as.list ...
}


run_model <- function(parameters, scenario, time, initcondrun) {
  # ************************************************************************************* 
  
  # all initial conditions must be integers
  initoderun <- initcondrun
  staterun <- initoderun
  inp <- inputs(parameters, scenario)
  
  # Solve the ODEs and get output
  timesrun <- time  # Model run time
  
  #Solve ODE
  outoderun <- ode(y=staterun, times=timesrun, func=epiModel, 
                   parms=parameters, 
                   input=inp, scenario=scenario)
  
  # Compute transitions at each time step
  tranoderun<-matrix(0,nrow=dim(outoderun)[1],ncol=N*A)
  #print(dim(tranoderun))
  # for (ti in 1:(tsteps+1)){
  #   #print(ti)
  #   tranoderun[ti,]<-t(malrates(outoderun[ti,2:(1+V)],inp,parameters,migrationMatrix, 0,ti,scenario ))
  # }
  #Compute outputs
  #ppout<-postproc(parameters,outoderun,tranoderun)
  
  modeltimes<-outoderun[,1]+startyear
  
  #inc_pred_ode<-ppout[,1:N]
  
  
  MALout<-list(#as.data.table(inc_pred_ode), #2
    outoderun #3
  )
  
  # state_names <- c()
  # for (n in 1:N) state_names <- c(state_names, paste0(var_names, "_p", n))
  # colnames(outoderun) <- c("time", state_names)
  # return(outoderun)
  
  return(MALout)
}



scenario <- NULL

initcondrun <- NULL
for (n in 1:N) {
  initcondrun <- c(initcondrun, c(rep(0,(B)))); 
  initcondrun[varind[1,n]] <- 0 # Immune Children
  initcondrun[varind[2,n]] <- 1000# Susceptible Childre
  initcondrun[varind[3,n]] <- 100#  Children Exposed Vivax
  initcondrun[varind[14,n]] <- 100; # Nothing in Recovered initially
  
  
  for(j in 5:85)
  {
    initcondrun[varind[j, n]] = 0
    
  }
  initcondrun[varind[86, n]] = 10000
  initcondrun[varind[87, n]] = 50000
  initcondrun[varind[88, n]] = 0
  initcondrun[varind[89, n]] = 50000
  initcondrun[varind[90, n]] = 0
}

initcondrun[varind[1,4]] <- 0 # Immune Children
initcondrun[varind[2,4]] <- 2500# Susceptible Childre
initcondrun[varind[3,4]] <- 2000#  Children Exposed Vivax
initcondrun[varind[14,4]] <- 50; # Nothing in Recovered initially

initcondrun[varind[1,11]] <- 0 # Immune Children
initcondrun[varind[2,11]] <- 50000# Susceptible Childre
initcondrun[varind[3,11]] <- 10000#  Children Exposed Vivax
initcondrun[varind[14,11]] <- 25000; # Nothing in Recovered initially

tryCatch({
  mo <- run_model(parameters, scenario, time, initcondrun)
}, error = function(e) {
  message("Model failed: ", conditionMessage(e))
  traceback()
})


col_idx <- function(vname, patch) {
  j <- match(vname, var_names)     # e.g., vname = "S_ch"
  if (is.na(j)) stop("Unknown variable name: ", vname)
  1 + varind[j, patch]             # +1 for the time column
}

# vector of that variable over time
get_series <- function(out, vname, patch) {
  out[, col_idx(vname, patch), drop = TRUE]
}


out <- mo[[1]]  # deSolve output matrix

get_series <- function(out, vname, patch) {
  out[, col_idx(vname, patch), drop = TRUE]
}

S_ch_p5       <- get_series(out, "S_ch", 5)
Iv_RDT_TP_p5  <- get_series(out, "Iv_RDT_TP_ch", 5)
Iv_RDT_FN_p5  <- get_series(out, "Iv_RDT_FN_ch", 5)
Iv_M_TP_p5    <- get_series(out, "Iv_M_TP_ch", 5)
Iv_M_FN_p5    <- get_series(out, "Iv_M_FN_ch", 5)
Av_p5         <- get_series(out, "Av_ch", 5)

IF_RDT_TP_p5  <- get_series(out, "IF_RDT_TP_ch", 5)
IF_RDT_FN_p5  <- get_series(out, "IF_RDT_FN_ch", 5)
IF_M_TP_p5    <- get_series(out, "IF_M_TP_ch", 5)
IF_M_FN_p5    <- get_series(out, "IF_M_FN_ch", 5)
AF_p5         <- get_series(out, "AF_ch", 5)

totalVivaxInfection      <- Iv_RDT_TP_p5 + Iv_RDT_FN_p5 + Iv_M_TP_p5 + Iv_M_FN_p5 + Av_p5
totalFalciparumInfection <- IF_RDT_TP_p5 + IF_RDT_FN_p5 + IF_M_TP_p5 + IF_M_FN_p5 + AF_p5

plot(out[,1], S_ch_p5, type = "l")   # use solver’s time column
lines(out[,1], totalVivaxInfection, col = "red")
lines(out[,1], totalFalciparumInfection, col = "blue")
