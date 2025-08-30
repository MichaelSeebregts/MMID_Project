library(deSolve)

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
  void EQ(NumericVector eq, NumericVector transit, 
          IntegerVector transitionsiu1, IntegerVector transitionsiu2, 
          IntegerVector transitionsiv1, IntegerVector transitionsiv2) {
    int i, iu1, iu2, iv1, iv2;
    // Zero out eq
    eq.fill(0.0);
    // Fill equations with new deltas
    //Rcpp::Rcout << transit.length() << std::endl;
    for(i=0; i<transit.length(); i++) {
      iu1 = transitionsiu1[i]-1;
      iv1 = transitionsiv1[i];
      iu2 = transitionsiu2[i]-1;
      iv2 = transitionsiv2[i];
      
      // Include transition differential in the relevant eq components
      
      eq[iu1] += transit[i]*iv1;
      eq[iu2] += transit[i]*iv2;
      
      //Rcpp::Rcout << "eq["<<iu1<<"] += "<<transit[i]<<"*"<<iv1<<";" << std::endl;
      //Rcpp::Rcout << "eq["<<iu2<<"] += "<<transit[i]<<"*"<<iv2<<";" << std::endl;
    }
  }
')


N = 12
B = 22 + 21*3 + 5  # number of variables per patch
A = 38  + # Transitions for Child
    36*3 + # Transitions for Men, Non Pregnant Women and Pregnant Women
    ((22 + 21*3)*(N-1)) + # Migration 
    (22 + 21*3) + # Deaths
    (21 + 21 + 21*2) + # Aging for children and women moving from not pregnant to pregnant 
    (6)*4 + # Death when infected with vivax or falciparum 
    5 + # Mosquito movements
    5 # Mosquito Deaths
V = N*B
L = N*A

startyear = 2022 # starting year of simulation 2022-01-01

tyears = 1 # total weeks of simulation 
dtout = 7/365 # output timestep
tsteps = round(tyears/dtout) # number of time steps
time = startyear+seq(0, tyears, dtout) # time vector

alivepopHum = c(1:(B - 5))
alivepopMos = c((B - 4):B)
popc <- rep(0,N)
popM <- rep(0,N)

childClass = c("ch")
otherClasses = c("m","wnp","wp")
compartmentsChildren = c("IT", "S", "Ev", "Av", "Iv_RDT_TP", "Iv_RDT_FN", "Iv_M_TP", "Iv_M_FN", "Tv", "notTv", "Tv_D", "Ev_D", "Rv", "EF", "AF", "IF_RDT_TP", "IF_RDT_FN", "IF_M_TP", "IF_M_FN", "TF", "notTF", "RF")
compartmentsOther = c("S", "Ev", "Av", "Iv_RDT_TP", "Iv_RDT_FN", "Iv_M_TP", "Iv_M_FN", "Tv", "notTv", "Tv_D", "Ev_D", "Rv", "EF", "AF", "IF_RDT_TP", "IF_RDT_FN", "IF_M_TP", "IF_M_FN", "TF", "notTF", "RF")
var_namesChildren = as.vector(sapply(childClass, function(cl) paste0(compartmentsChildren, "_", cl)))
var_namesOther = as.vector(sapply(otherClasses, function(cl) paste0(compartmentsOther, "_", cl)))

var_namesMosquitos = c("S_Mo", "Ev_Mo", "Iv_Mo", "EF_Mo", "IF_Mo")

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


transitions = matrix(0, nrow=L, ncol=4)
for (n in 1:N) {
  indAdd = 0
  
  migrationVector = c(1:12)[-n]
  
  # Birth
  transitions[traind[1 + indAdd,n],] <- c(varind[1,n],  0, varind[1,n], +1) #    -> IT
  indAdd = indAdd +1
  
  # Loss of immunity
  transitions[traind[2 + indAdd,n],] <- c(varind[1,n],  -1, varind[2,n], +1) #  IT  -> S
  indAdd = indAdd + 1
  
  print("immunity Done")
  
  # Death
  for (i in 1:B)
  {
    transitions[traind[indAdd + i,n],] <- c(varind[i,n],  -1, varind[1,n], 0) #  X -> 
  }
  indAdd = indAdd + B
  print(indAdd)
  
  print("death Done")
  
  # Vivax 
  
  for (j in 1:4)
  {
    
    if (j == 1)
    {
      transitions[traind[1 + indAdd, n],] <- c(varind[4,n], +1, varind[3,n], -1)
      transitions[traind[2 + indAdd, n],] <- c(varind[5,n], +1, varind[3,n], -1)
      transitions[traind[3 + indAdd, n],] <- c(varind[6,n], +1, varind[3,n], -1)
      transitions[traind[4 + indAdd, n],] <- c(varind[7,n], +1, varind[3,n], -1)
      transitions[traind[5 + indAdd, n],] <- c(varind[8,n], +1, varind[3,n], -1)
      transitions[traind[6 + indAdd, n],] <- c(varind[4,n], -1, varind[12,n], +1)
      transitions[traind[7 + indAdd, n],] <- c(varind[5,n], -1, varind[9,n], +1)
      transitions[traind[8 + indAdd, n],] <- c(varind[7,n], -1, varind[9,n], +1)
      transitions[traind[9 + indAdd, n],] <- c(varind[6,n], -1, varind[12,n], +1)
      transitions[traind[10 + indAdd, n],] <- c(varind[8,n], -1, varind[12,n], +1)
      transitions[traind[11 + indAdd, n],] <- c(varind[9,n], -1, varind[10,n], +1)
      transitions[traind[12 + indAdd, n],] <- c(varind[9,n], -1, varind[11,n], +1)
      transitions[traind[13 + indAdd, n],] <- c(varind[10,n], -1, varind[12,n], +1)
      transitions[traind[14 + indAdd, n],] <- c(varind[13,n], +1, varind[11,n], -1)
      transitions[traind[15 + indAdd, n],] <- c(varind[12,n], +1, varind[11,n], -1)
      transitions[traind[16 + indAdd, n],] <- c(varind[4,n], +1, varind[12,n], -1)
      transitions[traind[17 + indAdd, n],] <- c(varind[5,n], +1, varind[12,n], -1)
      transitions[traind[18 + indAdd, n],] <- c(varind[6,n], +1, varind[12,n], -1)
      transitions[traind[19 + indAdd, n],] <- c(varind[7,n], +1, varind[12,n], -1)
      transitions[traind[20 + indAdd, n],] <- c(varind[8,n], +1, varind[12,n], -1)
      transitions[traind[21 + indAdd, n],] <- c(varind[2,n], +1, varind[13,n], -1)
      transitions[traind[22 + indAdd, n],] <- c(varind[3,n], +1, varind[2,n], -1)
      
      indAdd = indAdd + 22
      print(indAdd)
      
    }
    
    if (j > 1)
    {
      transitions[traind[1 + indAdd, n],]  <- c(varind[3 + 22*(j-1), n],  +1, varind[2 + 22*(j-1), n],  -1)
      transitions[traind[2 + indAdd, n],]  <- c(varind[4 + 22*(j-1), n],  +1, varind[2 + 22*(j-1), n],  -1)
      transitions[traind[3 + indAdd, n],]  <- c(varind[5 + 22*(j-1), n],  +1, varind[2 + 22*(j-1), n],  -1)
      transitions[traind[4 + indAdd, n],]  <- c(varind[6 + 22*(j-1), n],  +1, varind[2 + 22*(j-1), n],  -1)
      transitions[traind[5 + indAdd, n],]  <- c(varind[7 + 22*(j-1), n],  +1, varind[2 + 22*(j-1), n],  -1)
      transitions[traind[6 + indAdd, n],]  <- c(varind[3 + 22*(j-1), n],  -1, varind[11 + 22*(j-1), n], +1)
      transitions[traind[7 + indAdd, n],]  <- c(varind[4 + 22*(j-1), n],  -1, varind[8 + 22*(j-1), n],  +1)
      transitions[traind[8 + indAdd, n],] <- c(varind[6 + 22*(j-1), n],  -1, varind[8 + 22*(j-1), n],  +1)
      transitions[traind[9 + indAdd, n],] <- c(varind[5 + 22*(j-1), n],  -1, varind[11 + 22*(j-1), n], +1)
      transitions[traind[10 + indAdd, n],] <- c(varind[7 + 22*(j-1), n],  -1, varind[11 + 22*(j-1), n], +1)
      transitions[traind[11 + indAdd, n],] <- c(varind[8 + 22*(j-1), n],  -1, varind[9 + 22*(j-1), n], +1)
      transitions[traind[12 + indAdd, n],] <- c(varind[8 + 22*(j-1), n],  -1, varind[10 + 22*(j-1), n], +1)
      transitions[traind[13 + indAdd, n],] <- c(varind[9 + 22*(j-1), n], -1, varind[11 + 22*(j-1), n], +1)
      transitions[traind[14 + indAdd, n],] <- c(varind[12 + 22*(j-1), n], +1, varind[10 + 22*(j-1), n], -1)
      transitions[traind[15 + indAdd, n],] <- c(varind[11 + 22*(j-1), n], +1, varind[10 + 22*(j-1), n], -1)
      transitions[traind[16 + indAdd, n],] <- c(varind[3 + 22*(j-1), n],  +1, varind[11 + 22*(j-1), n], -1)
      transitions[traind[17 + indAdd, n],] <- c(varind[4 + 22*(j-1), n],  +1, varind[11 + 22*(j-1), n], -1)
      transitions[traind[18 + indAdd, n],] <- c(varind[5 + 22*(j-1), n],  +1, varind[11 + 22*(j-1), n], -1)
      transitions[traind[19 + indAdd, n],] <- c(varind[6 + 22*(j-1), n],  +1, varind[11 + 22*(j-1), n], -1)
      transitions[traind[20 + indAdd, n],] <- c(varind[7 + 22*(j-1), n],  +1, varind[11 + 22*(j-1), n], -1)
      transitions[traind[21 + indAdd, n],] <- c(varind[1 + 22*(j-1), n],  +1, varind[12 + 22*(j-1), n], -1)
      transitions[traind[22 + indAdd, n],] <- c(varind[2 + 22*(j-1), n],  +1, varind[1 + 22*(j-1), n], -1)
      indAdd = indAdd + 22
      print(indAdd)
    }
    
  }
  print(indAdd)
  
  print("vivax Done")
  # Falciparum
  
  for (j in 1:4)
  {
    

    if (j == 1)
    {
      print(indAdd)
      transitions[traind[1 + indAdd, n],] <- c(varind[14, n],  +1, varind[2, n],  -1)
      transitions[traind[2 + indAdd, n],] <- c(varind[15, n],  +1, varind[14, n], -1)
      transitions[traind[3 + indAdd, n],] <- c(varind[16, n],  +1, varind[14, n], -1)
      transitions[traind[4 + indAdd, n],] <- c(varind[17, n],  +1, varind[14, n], -1)
      transitions[traind[5 + indAdd, n],] <- c(varind[18, n],  +1, varind[14, n], -1)
      transitions[traind[6 + indAdd, n],] <- c(varind[19, n],  +1, varind[14, n], -1)
      transitions[traind[7 + indAdd, n],] <- c(varind[20, n],  +1, varind[16, n], -1)
      transitions[traind[8 + indAdd, n],] <- c(varind[20, n],  +1, varind[18, n], -1)
      transitions[traind[9 + indAdd, n],] <- c(varind[22, n],  +1, varind[17, n], -1)
      transitions[traind[10 + indAdd, n],] <- c(varind[22, n],  +1, varind[19, n], -1)
      transitions[traind[11 + indAdd, n],] <- c(varind[22, n],  +1, varind[15, n], -1)
      transitions[traind[12 + indAdd, n],] <- c(varind[22, n],  +1, varind[21, n], -1)
      transitions[traind[13 + indAdd, n],] <- c(varind[22, n],  +1, varind[20, n], -1)
      transitions[traind[14 + indAdd, n],] <- c(varind[2, n],   +1, varind[22, n], -1)
      indAdd = indAdd + 14
      
      

    }
    
    if (j > 1)
    {
      
      transitions[traind[1 + indAdd, n],] <- c(varind[13+ 22 + 21*(j-2), n],  +1, varind[1+ 22 + 21*(j-2), n],  -1)
      transitions[traind[2 + indAdd, n],] <- c(varind[14+ 22 + 21*(j-2), n],  +1, varind[13+ 22 + 21*(j-2), n], -1)
      transitions[traind[3 + indAdd, n],] <- c(varind[15+ 22 + 21*(j-2), n],  +1, varind[13+ 22 + 21*(j-2), n], -1)
      transitions[traind[4 + indAdd, n],] <- c(varind[16+ 22 + 21*(j-2), n],  +1, varind[13+ 22 + 21*(j-2), n], -1)
      transitions[traind[5 + indAdd, n],] <- c(varind[17+ 22 + 21*(j-2), n],  +1, varind[13+ 22 + 21*(j-2), n], -1)
      transitions[traind[6 + indAdd, n],] <- c(varind[18+ 22 + 21*(j-2), n],  +1, varind[13+ 22 + 21*(j-2), n], -1)
      transitions[traind[7 + indAdd, n],] <- c(varind[19+ 22 + 21*(j-2), n],  +1, varind[15+ 22 + 21*(j-2), n], -1)
      transitions[traind[8 + indAdd, n],] <- c(varind[19+ 22 + 21*(j-2), n],  +1, varind[17+ 22 + 21*(j-2), n], -1)
      transitions[traind[9 + indAdd, n],] <- c(varind[21+ 22 + 21*(j-2), n],  +1, varind[16+ 22 + 21*(j-2), n], -1)
      transitions[traind[10 + indAdd, n],] <- c(varind[21+ 22 + 21*(j-2), n],  +1, varind[18+ 22 + 21*(j-2), n], -1)
      transitions[traind[11 + indAdd, n],] <- c(varind[21+ 22 + 21*(j-2), n],  +1, varind[14+ 22 + 21*(j-2), n], -1)
      transitions[traind[12 + indAdd, n],] <- c(varind[21+ 22 + 21*(j-2), n],  +1, varind[20+ 22 + 21*(j-2), n], -1)
      transitions[traind[13 + indAdd, n],] <- c(varind[21+ 22 + 21*(j-2), n],  +1, varind[19+ 22 + 21*(j-2), n], -1)
      transitions[traind[14 + indAdd, n],] <- c(varind[1+ 22 + 21*(j-2), n],   +1, varind[21+ 22 + 21*(j-2), n], -1)
      indAdd = indAdd + 14
    }
  }
  print("falciparum Done")
  print(indAdd)
  # Aging
  
  k = 1
  
  for (i in 23:43)
  {
    transitions[traind[k + indAdd,n],] <- c(varind[i-21,n],  -1, varind[1,n], +1) # Aging from Child to Man
    indAdd = indAdd + 1
    k = k + 1
  }
  
  print(indAdd)
  k = 1
  for (i in 44:64)
  {
    transitions[traind[k + indAdd,n],] <- c(varind[i-42,n],  -1, varind[1,n], +1) # Aging from Child to Woman 
    indAdd = indAdd + 1
    k = k + 1
  }
  print(indAdd)
  print("aging Done")
  
  # Pregnancy
  print(indAdd)
  k = 1
  for (i in 44:64)
  {
    transitions[traind[k + indAdd,n],] <- c(varind[i,n],  -1, varind[i + 21,n], +1) # Not Pregnant to Pregnant
    indAdd = indAdd + 1
    k = k + 1
  }
  
  print(indAdd)
  k = 1

  for (i in 65:85)
  {
    transitions[traind[k + indAdd,n],] <- c(varind[i,n],  -1, varind[i-21,n], +1) # Pregnant to Not Pregnant
    indAdd = indAdd + 1
    k = k + 1
  }
  
  print("pregnancy Done")
  
  # Migration
  print(indAdd)
  k = 1
  for (z in 1:11)
  {
    for (i in 1:(B-5))
    {
      transitions[traind[indAdd + 1,n],] <- c(varind[i,n],  -1, varind[i,migrationVector[z]], +1) #  Migrating from n to every 1 to 12, excluding n
      #print(indAdd + k)
      indAdd = indAdd + 1
    }
    
  }
  
  print("migration Done")
  
  # Deaths from Vivax and Falciparum
  print(indAdd)
  transitions[traind[1 + indAdd,n],] <- c(varind[4,n],  -1, varind[1,n], 0)
  transitions[traind[2 + indAdd,n],] <- c(varind[5,n],  -1, varind[1,n], 0)
  transitions[traind[3 + indAdd,n],] <- c(varind[7,n],  -1, varind[1,n], 0)
  transitions[traind[4 + indAdd,n],] <- c(varind[15,n],  -1, varind[1,n], 0)
  transitions[traind[5 + indAdd,n],] <- c(varind[16,n],  -1, varind[1,n], 0)
  transitions[traind[6 + indAdd,n],] <- c(varind[18,n],  -1, varind[1,n], 0)
  
  transitions[traind[7 + indAdd,n],] <- c(varind[26,n],  -1, varind[1,n], 0)
  transitions[traind[8 + indAdd,n],] <- c(varind[27,n],  -1, varind[1,n], 0)
  transitions[traind[9 + indAdd,n],] <- c(varind[29,n],  -1, varind[1,n], 0)
  transitions[traind[10 + indAdd,n],] <- c(varind[37,n],  -1, varind[1,n], 0)
  transitions[traind[11 + indAdd,n],] <- c(varind[38,n],  -1, varind[1,n], 0)
  transitions[traind[12 + indAdd,n],] <- c(varind[40,n],  -1, varind[1,n], 0)
  
  transitions[traind[13 + indAdd,n],] <- c(varind[47,n],  -1, varind[1,n], 0)
  transitions[traind[14 + indAdd,n],] <- c(varind[48,n],  -1, varind[1,n], 0)
  transitions[traind[15 + indAdd,n],] <- c(varind[50,n],  -1, varind[1,n], 0)
  transitions[traind[16 + indAdd,n],] <- c(varind[58,n],  -1, varind[1,n], 0)
  transitions[traind[17 + indAdd,n],] <- c(varind[59,n],  -1, varind[1,n], 0)
  transitions[traind[18 + indAdd,n],] <- c(varind[61,n],  -1, varind[1,n], 0)
  
  transitions[traind[19 + indAdd,n],] <- c(varind[68,n],  -1, varind[1,n], 0)
  transitions[traind[20 + indAdd,n],] <- c(varind[69,n],  -1, varind[1,n], 0)
  transitions[traind[21 + indAdd,n],] <- c(varind[71,n],  -1, varind[1,n], 0)
  transitions[traind[22 + indAdd,n],] <- c(varind[79,n],  -1, varind[1,n], 0)
  transitions[traind[23 + indAdd,n],] <- c(varind[80,n],  -1, varind[1,n], 0)
  transitions[traind[24 + indAdd,n],] <- c(varind[82,n],  -1, varind[1,n], 0)
  
  print("Death from Vivax and Falciparum Done")
  
  # Mosquito
  
  transitions[traind[25 + indAdd,n],] <- c(varind[86,n],  0, varind[86,n], 1) # Birth
  # Mosquito deaths already included above
  
  transitions[traind[26 + indAdd,n],] <- c(varind[86,n],  -1, varind[87,n], 1) # S -> Ev
  transitions[traind[27 + indAdd,n],] <- c(varind[86,n],  -1, varind[89,n], 1) # S -> Ef
  
  transitions[traind[28 + indAdd,n],] <- c(varind[87,n],  -1, varind[88,n], 1) # Ev -> Iv
  transitions[traind[29 + indAdd,n],] <- c(varind[89,n],  -1, varind[90,n], 1) # Ef -> If
  
  print(indAdd + 29)
  
  print("Mosquitos Done")
    
  
}

transitionsiu1 <- transitions[,1]
transitionsiv1 <- transitions[,2]
transitionsiu2 <- transitions[,3]
transitionsiv2 <- transitions[,4]

inputs<-function(parameters, scenario){
  
}

parsHuman = readxl::read_excel("Parameters.xlsx", sheet = "Humans") %>% 
  pull(var=Value, name=Name)
migrationMatrix = readxl::read_excel("Parameters.xlsx", sheet = "migrationMatrix")
parsMosquito = readxl::read_excel("Parameters.xlsx", sheet = "Mosquitos") %>% 
  pull(var=Value, name=Name)

parameters = c(parsHuman, parsMosquito)



get_comp_vecs <- function(names, suffix, env = parent.frame()) {
  sapply(names, function(nm) rlang::eval_tidy(rlang::sym(paste0(nm, "_", suffix)), env = env))
}

malrates <- function(x, input, parameters, matrix, t, ti, scenario) {
  with(as.list(c(parameters, scenario, getStates(x))), {
    
    lambda_V = 0.8
    lambda_F = 0.8
    lambda_v_ch = 0
    lambda_F_ch = 0
    lambda_v_wp = 0 
    lambda_F_wp = 0 
    notPregnant = 1/40
    m1 = 0
    m2 = 0
    m3 = 0
    m4 = 0
    m5 = 0
    m6 = 0
    m7 = 0
    m8 = 0
    m9 = 0
    m10 = 0
    m11 = 0
    
    m1_m = 0
    m2_m = 0
    m3_m = 0
    m4_m = 0
    m5_m = 0
    m6_m = 0
    m7_m = 0
    m8_m = 0
    m9_m = 0
    m10_m = 0
    m11_m = 0
    
    migrationFor_ch_wp_wnp = c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11)
    migrationFor_m = c(m1_m, m2_m, m3_m, m4_m, m5_m, m6_m, m7_m, m8_m, m9_m, m10_m, m11_m)
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print(S_ch)
  
    
    emptyCat_ch = c(IT_ch, S_ch, Ev_ch, Av_ch, Iv_RDT_TP_ch, Iv_RDT_FN_ch,
                 Iv_M_TP_ch, Iv_M_FN_ch, Tv_ch, notTv_ch, Tv_D_ch, Ev_D_ch, Rv_ch, 
                 EF_ch, AF_ch, IF_RDT_TP_ch, IF_RDT_FN_ch, IF_M_TP_ch, 
                 IF_M_FN_ch, TF_ch, notTF_ch, RF_ch)
    emptyCat_m = c(S_m, Ev_m, Av_m, Iv_RDT_TP_m, Iv_RDT_FN_m,
                   Iv_M_TP_m, Iv_M_FN_m, Tv_m, notTv_m, Tv_D_m, Ev_D_m, Rv_m, 
                   EF_m, AF_m, IF_RDT_TP_m, IF_RDT_FN_m, IF_M_TP_m, 
                   IF_M_FN_m, TF_m, notTF_m, RF_m)
    emptyCat_wp = c( S_wp, Ev_wp, Av_wp, Iv_RDT_TP_wp, Iv_RDT_FN_wp,
                   Iv_M_TP_wp, Iv_M_FN_wp, Tv_wp, notTv_wp, Tv_D_wp, Ev_D_wp, Rv_wp, 
                   EF_wp, AF_wp, IF_RDT_TP_wp, IF_RDT_FN_wp, IF_M_TP_wp, 
                   IF_M_FN_wp, TF_wp, notTF_wp, RF_wp)
    emptyCat_wnp = c(S_wnp, Ev_wnp, Av_wnp, Iv_RDT_TP_wnp, Iv_RDT_FN_wnp,
                   Iv_M_TP_wnp, Iv_M_FN_wnp, Tv_wnp, notTv_wnp, Tv_D_wnp, Ev_D_wnp, Rv_wnp, 
                   EF_wnp, AF_wnp, IF_RDT_TP_wnp, IF_RDT_FN_wnp, IF_M_TP_wnp, 
                   IF_M_FN_wnp, TF_wnp, notTF_wnp, RF_wnp)
    
    #emptyCat = c(S, Ev, Av, Iv_RDT_TP, Iv_RDT_FN,
     #                Iv_M_TP, Iv_M_FN, Tv, notTv, Tv_D, Ev_D, Rv, 
    #                 EF, AF, IF_RDT_TP, IF_RDT_FN, IF_M_TP, 
      #               IF_M_FN, TF, notTF, RF)
    
# compartmentsChildren = c("IT", "S", "Ev", "Av", "Iv_RDT_TP", "Iv_RDT_FN", "Iv_M_TP", "Iv_M_FN", 
 #   "Tv", "notTv", "Tv_D", "Ev_D", "Rv", "EF", "AF", "IF_RDT_TP", "IF_RDT_FN", "IF_M_TP", 
 #   "IF_M_FN", "TF", "notTF", "RF")    
    
  # var_namesMosquitos = c("S_Mo", "Ev_Mo", "Iv_Mo", "EF_Mo", "IF_Mo")
    birth = c(b*.pop)
    loss_immune = c(eta*IT_ch)
    
    # death_ch = c(mu*IT_ch, mu*S_ch, mu*Ev_ch, mu_ch_V*Av_ch, mu_ch_V*Iv_RDT_TP_ch, mu_ch_V*Iv_RDT_FN_ch,
    #              mu_ch_V*Iv_M_TP_ch, mu_ch_V*Iv_M_FN_ch, mu*Tv_ch, mu*notTv_ch, mu*Tv_D_ch, mu*Ev_D_ch, mu*Rv_ch, 
    #              mu*EF_ch, mu_ch_F*AF_ch, mu_ch_F*IF_RDT_TP_ch, mu_ch_F*IF_RDT_FN_ch, mu_ch_F*IF_M_TP_ch, 
    #              mu_ch_F*IF_M_FN_ch, mu*TF_ch, mu*notTF_ch, mu*RF_ch)
    # death_m_wnp_cat = c(mu*IT, mu*S, mu*Ev, mu_V*Av, mu_V*Iv_RDT_TP, mu_V*Iv_RDT_FN,
    #                 mu_V*Iv_M_TP, mu_V*Iv_M_FN, mu*Tv, mu*notTv, mu*Tv_D, mu*Ev_D, mu*Rv, 
    #                 mu*EF, mu_F*AF, mu_F*IF_RDT_TP, mu_F*IF_RDT_FN, mu_F*IF_M_TP, 
    #                 mu_F*IF_M_FN, mu*TF, mu*notTF, mu*RF)
    
    death_ch = mu*emptyCat_ch
    death_m = mu*emptyCat_m
    death_wp = mu*emptyCat_wp
    death_wnp = mu*emptyCat_wnp
    
    
    # death_m = suffix_exprs(death_m_wnp_cat, "m")
    # death_wnp = suffix_exprs(death_m_wnp_cat, "wnp")
    # death_m_wnp = asvector(sapply(c("m", "wp"), function(cl) paste0(death_m_wnp, "_", cl)))
    # 
    # death_wp = c(mu*IT_wp, mu*S_wp, mu*Ev_wp, mu_wp_V*Av_wp, mu_wp_V*Iv_RDT_TP_wp, mu_wp_V*Iv_RDT_FN_wp,
    #              mu_wp_V*Iv_M_TP_wp, mu_wp_V*Iv_M_FN_wp, mu*Tv_wp, mu*notTv_wp, mu*Tv_D_wp, mu*Ev_D_wp, mu*Rv_wp, 
    #              mu*EF_wp, mu_wp_F*AF_wp, mu_wp_F*IF_RDT_TP_wp, mu_wp_F*IF_RDT_FN_wp, mu_wp_F*IF_M_TP_wp, 
    #              mu_wp_F*IF_M_FN_wp, mu*TF_wp, mu*notTF_wp, mu*RF_wp)
    
    deathMos = c(mu_mos*S_Mo, mu_mos*Ev_Mo, mu_mos*Iv_Mo, mu_mos*EF_Mo, mu_mos*IF_Mo)
    
    vivaxPath_ch = c(sigma_v*pA_v*Ev_ch, sigma_v*(1 - pA_v)*pTP_RDT*Ev_ch, sigma_v*(1 - pA_v)*pFN_RDT*Ev_ch, sigma_v*(1 - pA_v)*pTP_M*Ev_ch, 
                  sigma_v*(1 - pA_v)*pFN_RDT*Ev_ch, delta_v*Av_ch, tau*Iv_RDT_TP_ch, tau*Iv_M_TP_ch, delta_v*Iv_RDT_FN_ch, delta_v*Iv_M_FN_ch, 
                  (1-pCT)*r*Tv_ch, pCT*r*Tv_ch, phi*notTv_ch, (1-pTD)*theta*Tv_D_ch, pTD*theta*Tv_D_ch, smallOmega*pA_v*Ev_ch, 
                  smallOmega*(1 - pA_v)*pTP_RDT*Ev_ch, smallOmega*(1 - pA_v)*pFN_RDT*Ev_ch, smallOmega*(1 - pA_v)*pTP_M*Ev_ch, 
                  smallOmega*(1 - pA_v)*pFN_RDT*Ev_ch, rho*Rv_ch, lambda_v*S_ch)
    
    vivaxPath_m = c(sigma_v*pA_v*Ev_m, sigma_v*(1 - pA_v)*pTP_RDT*Ev_m, sigma_v*(1 - pA_v)*pFN_RDT*Ev_m, sigma_v*(1 - pA_v)*pTP_M*Ev_m, 
                     sigma_v*(1 - pA_v)*pFN_RDT*Ev_m, delta_v*Av_m, tau*Iv_RDT_TP_m, tau*Iv_M_TP_m, delta_v*Iv_RDT_FN_m, delta_v*Iv_M_FN_m, 
                     (1-pCT)*r*Tv_m, pCT*r*Tv_m, phi*notTv_m, (1-pTD)*theta*Tv_D_m, pTD*theta*Tv_D_m, smallOmega*pA_v*Ev_m, 
                     smallOmega*(1 - pA_v)*pTP_RDT*Ev_m, smallOmega*(1 - pA_v)*pFN_RDT*Ev_m, smallOmega*(1 - pA_v)*pTP_M*Ev_m, 
                     smallOmega*(1 - pA_v)*pFN_RDT*Ev_m, rho*Rv_m, lambda_v*S_m)
    
    vivaxPath_wp = c(sigma_v*pA_v*Ev_wp, sigma_v*(1 - pA_v)*pTP_RDT*Ev_wp, sigma_v*(1 - pA_v)*pFN_RDT*Ev_wp, sigma_v*(1 - pA_v)*pTP_M*Ev_wp, 
                     sigma_v*(1 - pA_v)*pFN_RDT*Ev_wp, delta_v*Av_wp, tau*Iv_RDT_TP_wp, tau*Iv_M_TP_wp, delta_v*Iv_RDT_FN_wp, delta_v*Iv_M_FN_wp, 
                     (1-pCT)*r*Tv_wp, pCT*r*Tv_wp, phi*notTv_wp, (1-pTD)*theta*Tv_D_wp, pTD*theta*Tv_D_wp, smallOmega*pA_v*Ev_wp, 
                     smallOmega*(1 - pA_v)*pTP_RDT*Ev_wp, smallOmega*(1 - pA_v)*pFN_RDT*Ev_wp, smallOmega*(1 - pA_v)*pTP_M*Ev_wp, 
                     smallOmega*(1 - pA_v)*pFN_RDT*Ev_wp, rho*Rv_wp, lambda_v*S_wp)
    
    vivaxPath_wnp = c(sigma_v*pA_v*Ev_wnp, sigma_v*(1 - pA_v)*pTP_RDT*Ev_wnp, sigma_v*(1 - pA_v)*pFN_RDT*Ev_wnp, sigma_v*(1 - pA_v)*pTP_M*Ev_wnp, 
                     sigma_v*(1 - pA_v)*pFN_RDT*Ev_wnp, delta_v*Av_wnp, tau*Iv_RDT_TP_wnp, tau*Iv_M_TP_wnp, delta_v*Iv_RDT_FN_wnp, delta_v*Iv_M_FN_wnp, 
                     (1-pCT)*r*Tv_wnp, pCT*r*Tv_wnp, phi*notTv_wnp, (1-pTD)*theta*Tv_D_wnp, pTD*theta*Tv_D_wnp, smallOmega*pA_v*Ev_wnp, 
                     smallOmega*(1 - pA_v)*pTP_RDT*Ev_wnp, smallOmega*(1 - pA_v)*pFN_RDT*Ev_wnp, smallOmega*(1 - pA_v)*pTP_M*Ev_wnp, 
                     smallOmega*(1 - pA_v)*pFN_RDT*Ev_wnp, rho*Rv_wnp, lambda_v*S_wnp)
    
    
    
    falciparumPath_ch = c(lambda_F*S_ch, sigma_F*(1 - pA_F)*pTP_RDT*EF_ch, sigma_F*(1 - pA_F)*pFN_RDT*EF_ch, sigma_F*(1 - pA_F)*pTP_M*EF_ch, 
                       sigma_F*(1 - pA_F)*pFN_RDT*EF_ch, tau*IF_RDT_TP_ch, tau*IF_M_TP_ch, delta_F*IF_RDT_FN_ch, delta_F*IF_M_FN_ch, delta_F*AF_ch, 
                       delta_F*notTF_ch, r*TF_ch, rho*RF_ch) 
    falciparumPath_m = c(lambda_F*S_m, sigma_F*(1 - pA_F)*pTP_RDT*EF_m, sigma_F*(1 - pA_F)*pFN_RDT*EF_m, sigma_F*(1 - pA_F)*pTP_M*EF_m, 
                          sigma_F*(1 - pA_F)*pFN_RDT*EF_m, tau*IF_RDT_TP_m, tau*IF_M_TP_m, delta_F*IF_RDT_FN_m, delta_F*IF_M_FN_m, delta_F*AF_m, 
                          delta_F*notTF_m, r*TF_m, rho*RF_m)
    falciparumPath_wp = c(lambda_F*S_wp, sigma_F*(1 - pA_F)*pTP_RDT*EF_wp, sigma_F*(1 - pA_F)*pFN_RDT*EF_wp, sigma_F*(1 - pA_F)*pTP_M*EF_wp, 
                          sigma_F*(1 - pA_F)*pFN_RDT*EF_wp, tau*IF_RDT_TP_wp, tau*IF_M_TP_wp, delta_F*IF_RDT_FN_wp, delta_F*IF_M_FN_wp, delta_F*AF_wp, 
                          delta_F*notTF_wp, r*TF_wp, rho*RF_wp)
    falciparumPath_wnp = c(lambda_F*S_wnp, sigma_F*(1 - pA_F)*pTP_RDT*EF_wnp, sigma_F*(1 - pA_F)*pFN_RDT*EF_wnp, sigma_F*(1 - pA_F)*pTP_M*EF_wnp, 
                          sigma_F*(1 - pA_F)*pFN_RDT*EF_wnp, tau*IF_RDT_TP_wnp, tau*IF_M_TP_wnp, delta_F*IF_RDT_FN_wnp, delta_F*IF_M_FN_wnp, delta_F*AF_wnp, 
                          delta_F*notTF_wnp, r*TF_wnp, rho*RF_wnp)
    
    
    
    agingToMen = c(a*pMal*S_ch, a*pMal*Ev_ch, a*pMal*Av_ch, a*pMal*Iv_RDT_TP_ch, a*pMal*Iv_RDT_FN_ch,
                 a*pMal*Iv_M_TP_ch, a*pMal*Iv_M_FN_ch, a*pMal*Tv_ch, a*pMal*notTv_ch, a*pMal*Tv_D_ch, a*pMal*Ev_D_ch, a*pMal*Rv_ch, 
                 a*pMal*EF_ch, a*pMal*AF_ch, a*pMal*IF_RDT_TP_ch, a*pMal*IF_RDT_FN_ch, a*pMal*IF_M_TP_ch, 
                 a*pMal*IF_M_FN_ch, a*pMal*TF_ch, a*pMal*notTF_ch, a*pMal*RF_ch)
    agingToWomen = c(a*(1-pMal)*S_ch, a*(1-pMal)*Ev_ch, a*(1-pMal)*Av_ch, a*(1-pMal)*Iv_RDT_TP_ch, a*(1-pMal)*Iv_RDT_FN_ch,
                   a*(1-pMal)*Iv_M_TP_ch, a*(1-pMal)*Iv_M_FN_ch, a*(1-pMal)*Tv_ch, a*(1-pMal)*notTv_ch, a*(1-pMal)*Tv_D_ch, a*(1-pMal)*Ev_D_ch, a*(1-pMal)*Rv_ch, 
                   a*(1-pMal)*EF_ch, a*(1-pMal)*AF_ch, a*(1-pMal)*IF_RDT_TP_ch, a*(1-pMal)*IF_RDT_FN_ch, a*(1-pMal)*IF_M_TP_ch, 
                   a*(1-pMal)*IF_M_FN_ch, a*(1-pMal)*TF_ch, a*(1-pMal)*notTF_ch, a*(1-pMal)*RF_ch)
    # compartmentsChildren = c("IT", "S", "Ev", "Av", "Iv_RDT_TP", "Iv_RDT_FN", "Iv_M_TP", "Iv_M_FN", 
    #   "Tv", "notTv", "Tv_D", "Ev_D", "Rv", "EF", "AF", "IF_RDT_TP", "IF_RDT_FN", "IF_M_TP", 
    #   "IF_M_FN", "TF", "notTF", "RF") 
    
    wnpTowp = pregnant*emptyCat_wnp
    wpTownp = notPregnant*emptyCat_wp
    
    migrationVec = c()
    
    for (i in 1:11)
    {
      migrate_ch = migrationFor_ch_wp_wnp[i]*emptyCat_ch
      migrate_m = migrationFor_m[i]*emptyCat_m
      migrate_wp = migrationFor_ch_wp_wnp[i]*emptyCat_wp
      migrate_wnp = migrationFor_ch_wp_wnp[i]*emptyCat_wnp
      
      migrationVec = c(migrationVec, migrate_ch, migrate_m, migrate_wp, migrate_wnp)
    }
    
    additionalMortCompViv_ch = c(Av_ch, Iv_RDT_TP_ch, Iv_M_TP_ch)
    additionalMortCompFal_ch = c(AF_ch, IF_RDT_TP_ch, IF_M_TP_ch)
    
    additionalMortCompViv_m = c(Av_m, Iv_RDT_TP_m, Iv_M_TP_m)
    additionalMortCompFal_m = c(AF_m, IF_RDT_TP_m, IF_M_TP_m)
    
    additionalMortCompViv_wp = c(Av_wp, Iv_RDT_TP_wp, Iv_M_TP_wp)
    additionalMortCompFal_wp = c(AF_wp, IF_RDT_TP_wp, IF_M_TP_wp)
    
    additionalMortCompViv_wnp = c(Av_wnp, Iv_RDT_TP_wnp, Iv_M_TP_wnp)
    additionalMortCompFal_wnp = c(AF_wnp, IF_RDT_TP_wnp, IF_M_TP_wnp)
    
    addMortViv_ch = mu_ch_V*additionalMortCompViv_ch
    addMortFal_ch = mu_ch_F*additionalMortCompFal_ch
    
    addMortViv_m = mu_V*additionalMortCompViv_m
    addMortFal_m = mu_F*additionalMortCompFal_m
    
    addMortViv_wp = mu_wp_V*additionalMortCompViv_wp
    addMortFal_wp = mu_wp_F*additionalMortCompFal_wp
    
    addMortViv_wnp = mu_V*additionalMortCompViv_wnp
    addMortFal_wnp = mu_F*additionalMortCompFal_wnp
    
    # var_namesMosquitos = c("S_Mo", "Ev_Mo", "Iv_Mo", "EF_Mo", "IF_Mo")
    # Mosquito
    
    mosMovement = c(b_mos*.popMos, lambda_mos_V*S_Mo, gamma_mos_V*Ev_Mo, 
                    lambda_mos_F*S_Mo, gamma_mos_F*EF_Mo)
    
    tranrate <- array(c(
      birth, loss_immune, death_ch, death_m, death_wp, death_wnp, 
      deathMos, vivaxPath_ch, vivaxPath_m, vivaxPath_wp, vivaxPath_wnp, 
      falciparumPath_ch, falciparumPath_m, falciparumPath_wp, falciparumPath_wnp, 
      agingToMen, agingToWomen, wnpTowp, wpTownp, migrationVec, addMortViv_ch,
      addMortFal_ch, addMortViv_m, addMortFal_m, addMortViv_wp, addMortFal_wp, 
      addMortViv_wnp, addMortFal_wnp, mosMovement
      
      
      

      ), dim=c(N, A))
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


dZ <- rep(0.0, V) %>% as.numeric()
epiModel <- function(t, state, parameters, input, scenario) {
  with(as.list(c(state, parameters)), {
    # ************************************************************************************* #
    # define variables
    # ************************************************************************************* #
    
    Z <- state
    
    # rates of change
    ti <- 1
    transit <- malrates(state, input, parameters, migrationMatrix, t, ti, scenario)
    EQ(dZ, transit, transitionsiu1, transitionsiu2, transitionsiv1, transitionsiv2)
    # return the rate of change
    # print(c(t, dZ))
    # browser()
    list(c(dZ))
  }
  ) 
  # end with(as.list ...
}


run_model <- function(parameters, scenario) {
  # ************************************************************************************* 
  initcondrun <- NULL
  for (n in 1:N) {
    initcondrun <- c(initcondrun, c(rep(0,(B)))); 
    initcondrun[varind[1,n]] <- 0 # susceptible
    initcondrun[varind[2,n]] <- 2100# exposed
    initcondrun[varind[3,n]] <- 2900#  infectious
    initcondrun[varind[4,n]] <- 2900; # Nothing in Recovered initially
    for(j in 5:85)
    {
      initcondrun[varind[j, n]] = 0
      
    }
    initcondrun[varind[86, n]] = 10000
    initcondrun[varind[87, n]] = 50000
    initcondrun[varind[88, n]] = 50000
    initcondrun[varind[89, n]] = 0
    initcondrun[varind[90, n]] = 0
  }
  
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
  print("before tranoderun")
  tranoderun<-matrix(0,nrow=length(outoderun[,1]),ncol=length(transitions))
  for (ti in 1:(tsteps+1)){
    tranoderun[ti,]<-t(malrates(outoderun[ti,2:(1+V)],inp,parameters,migrationMatrix, 0,ti,scenario ))
  }
  #Compute outputs
  ppout<-postproc(parameters,outoderun,tranoderun)
  modeltimes<-outoderun[,1]+startyear
  
  inc_pred_ode<-ppout[,1:N]
  
  
  MALout<-list(as.data.table(inc_pred_ode), #2
               outoderun #3
               
  )
  
  return(MALout)
}



scenario <- NULL

tryCatch({
  mo <- run_model(parameters, scenario)
}, error = function(e) {
  message("Model failed: ", conditionMessage(e))
  traceback()
})

