
library(sampling)
library(survey)
library(data.table)
library(nonprobsvy)
library(jointCalib)
library(doSNOW)
library(progress)
library(foreach)
library(ggplot2)

source("codes/sim-code-yang-2020.R")
#save_file <- "results/paper-nonrob-results-20k-revised.rds"
#save_file <- "results/paper-nonrob-results-50k-revised.rds"
save_file <- "results/paper-nonrob-results-100k-revised.rds"

seed_number <- 2023-12-10
set.seed(seed_number)
N <- 100000
n <- 1000
x1 <- rnorm(N,1,1)
x2 <- rexp(N,1)
alp <- rnorm(N)
epsilon <- rnorm(N)
y11 <- 1 + x1 + x2 + alp + epsilon
y12 <- 0.5*(x1-1.5)^2 + x2^2 + alp + epsilon
y21 <- rbinom(N,1,plogis(1 + x1 + x2 + alp))
y22 <- rbinom(N,1,plogis(0.5*(x1-1.5)^2 + x2^2 + alp))
p1 <- plogis(x2)
p2 <- plogis(-3+(x1-1.5)^2+(x2-2)^2)
pop_data <- data.frame(x1,x2,y11,y12,y21,y22,p1,p2) |> setDT()
p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles2 <- seq(0.10, 0.90, 0.10)

formula_xs_p1a <- ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75
formula_xs_p1b <- ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75

formula_xs_p2a <- ~ x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +  x1_0.7 + x1_0.8 + 
  x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +  x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9

formula_xs_p2b <- ~ x1 + x2 + 
  x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +  x1_0.7 + x1_0.8 + x1_0.9 + 
  x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 + x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9

formula_y11_p1a <- as.formula(paste("y11", paste(as.character(formula_xs_p1a), collapse ="")))
formula_y11_p1b <- as.formula(paste("y11", paste(as.character(formula_xs_p1b), collapse ="")))
formula_y11_p2a <- as.formula(paste("y11", paste(as.character(formula_xs_p2a), collapse ="")))
formula_y11_p2b <- as.formula(paste("y11", paste(as.character(formula_xs_p2b), collapse ="")))

formula_y12_p1a <- as.formula(paste("y12", paste(as.character(formula_xs_p1a), collapse ="")))
formula_y12_p1b <- as.formula(paste("y12", paste(as.character(formula_xs_p1b), collapse ="")))
formula_y12_p2a <- as.formula(paste("y12", paste(as.character(formula_xs_p2a), collapse ="")))
formula_y12_p2b <- as.formula(paste("y12", paste(as.character(formula_xs_p2b), collapse ="")))

formula_y21_p1a <- as.formula(paste("y21", paste(as.character(formula_xs_p1a), collapse ="")))
formula_y21_p1b <- as.formula(paste("y21", paste(as.character(formula_xs_p1b), collapse ="")))
formula_y21_p2a <- as.formula(paste("y21", paste(as.character(formula_xs_p2a), collapse ="")))
formula_y21_p2b <- as.formula(paste("y21", paste(as.character(formula_xs_p2b), collapse ="")))

formula_y22_p1a <- as.formula(paste("y22", paste(as.character(formula_xs_p1a), collapse ="")))
formula_y22_p1b <- as.formula(paste("y22", paste(as.character(formula_xs_p1b), collapse ="")))
formula_y22_p2a <- as.formula(paste("y22", paste(as.character(formula_xs_p2a), collapse ="")))
formula_y22_p2b <- as.formula(paste("y22", paste(as.character(formula_xs_p2b), collapse ="")))

a <- Sys.time()
sims <- 600
cores <- 8
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation1 <- foreach(k=1:sims, 
                               .combine = rbind,
                               .packages = c("survey", "nonprobsvy", "jointCalib", "sampling", "data.table"),
                               .options.snow = opts,
                               .errorhandling = "remove") %dopar% {
                        yang_sim(k)                              
}
stopCluster(cl)
print(Sys.time() - a) 

setDT(results_simulation1)
results_simulation1[, c:=sum(is.nan(mean)), k]
results_simulation1[, d:=sum(lower_bound < 0, na.rm=T), k]
results_simulation1[y == "y11", true:=mean(y11)]
results_simulation1[y == "y12", true:=mean(y12)]
results_simulation1[y == "y21", true:=mean(y21)]
results_simulation1[y == "y22", true:=mean(y22)]

saveRDS(results_simulation1, save_file)
