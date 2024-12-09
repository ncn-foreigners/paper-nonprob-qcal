## plot 1

library(jointCalib)
library(latex2exp)
library(Cairo)
library(glue)
library(cobalt)

set.seed(20240109)
N <- 1000
x <- runif(N, 0, 80)
y <- 1300-(x-40)^2 + rnorm(N, 0, 300)
#y <- 1000+10*x+rnorm(N, 0, 300)
#y <- exp(-0.1 + 0.1*x) + rnorm(N, 0, 300)
probs2 <- seq(0.1, 0.9, 0.1)
quants_known1 <- list(x=quantile(x, c(0.25,0.5,0.75)))
quants_known2 <- list(x=quantile(x,probs2))
totals_known <- c(x=sum(x))
df <- data.frame(x, y)
x1 <- joint_calib_create_matrix(df, N, quants_known1)
colnames(x1) <- paste0("x_d", c(0.25, 0.50, 0.75))
colnames(x1) <- gsub("\\.", "", colnames(x1))
x2 <- joint_calib_create_matrix(df, N, quants_known2)
colnames(x2) <- paste0("x_p", probs2)
colnames(x2) <- gsub("\\.", "", colnames(x2))
df <- cbind(df, x1, x2)
m1 <- lm(y ~ x, df)
m2a <- lm(y ~ x_d025 + x_d05 + x_d075, df)
m2b <- lm(y ~ x + x_d025 + x_d05 + x_d075, df)
m3a <- lm(y ~ x_p01 + x_p02 + x_p03 + x_p04 + x_p05 + x_p06 + x_p07 + x_p08 + x_p09, df)
m3b <- lm(y ~ x + x_p01 + x_p02 + x_p03 + x_p04 + x_p05 + x_p06 + x_p07 + x_p08 + x_p09, df)
rsquares <- sapply(list(m1,m2a,m2b,m3a,m3b), function(x) round(summary(x)$r.square*100))
rsquares

cairo_pdf(file = "figs/fig1-example1.pdf", width = 13, height =5)
par(mfrow=c(1,3))
with(df, plot(x, y, main = TeX(r"($E(y_k|x_k)$,$R^2 \approx 0\%$)"), xlab = "X", ylab = "Y1", cex.main = 1.5))
with(df, points(x, fitted(m1), col = "red"))

#with(df, plot(x, y, main = TeX(r"($E(y_k|x_k)$,$R^2 \approx 49\%$)"), xlab = "X (quartiles used)", ylab = ""))
#with(df, points(x, fitted(m2a), col = "red"))

with(df, plot(x, y, 
              main = TeX(r"($E(y_k|x_k, a_k)$,$R^2 \approx 49\%$)"), xlab = "X (quartiles used)", ylab = "", 
              cex.main = 1.5))
with(df, points(x, fitted(m2b), col = "red"))

#with(df, plot(x, y, main = TeX(r"($E(y_k|x_k)$,$R^2 \approx 68\%$)"), xlab = "X (deciles used)", ylab = ""))
#with(df, points(x, fitted(m3a), col = "red"))

with(df, plot(x, y, main = TeX(r"($E(y_k|x_k, a_k)$,$R^2 \approx 68\%$)"), xlab = "X (deciles used)", 
              ylab = "", cex.main = 1.5))
with(df, points(x, fitted(m3b), col = "red"))

dev.off()

## fig2 -

set.seed(20240109)
N <- 1000
x <- rnorm(N,1,1)
p2 <- plogis(-3+(x-1.5)^2 + rnorm(N, 0, 0.5))
#p2 <- plogis(-3+x + rnorm(N, 0, 0.5))
y <- rbinom(N, 1, p2)
probs2 <- seq(0.1, 0.9, 0.1)
plot(x, p2)
quants_known1 <- list(x=quantile(x, c(0.25,0.5,0.75)))
quants_known2 <- list(x=quantile(x, probs2))
totals_known <- c(x=sum(x))
df <- data.frame(x, y)
x1 <- joint_calib_create_matrix(df, N, quants_known1)
colnames(x1) <- paste0("x_d", c(0.25, 0.50, 0.75))
colnames(x1) <- gsub("\\.", "", colnames(x1))
x2 <- joint_calib_create_matrix(df, N, quants_known2)
colnames(x2) <- paste0("x_p", probs2)
colnames(x2) <- gsub("\\.", "", colnames(x2))
df <- cbind(df, x1, x2)
m1 <- glm(y ~ x, df, family = binomial(link = "logit"))
m2a <- glm(y ~ x_d025 + x_d05 + x_d075, df, family = binomial(link = "logit"))
m2b <- glm(y ~ x + x_d025 + x_d05 + x_d075, df, family = binomial(link = "logit"))
m3a <- glm(y ~ x_p01 + x_p02 + x_p03 + x_p04 + x_p05 + x_p06 + x_p07 + x_p08 + x_p09, df, family = binomial(link = "logit"))
m3b <- glm(y ~ x + x_p01 + x_p02 + x_p03 + x_p04 + x_p05 + x_p06 + x_p07 + x_p08 + x_p09, df, family = binomial(link = "logit"))


cairo_pdf(file = "figs/fig1-example2.pdf", width = 13, height =5)
par(mfrow=c(1,3))
with(df, plot(x, p2, main = TeX(r"($P(y_k|x_k)$)"), xlab = "X", ylab = "P(Y2=1|X)", cex.main = 1.5))
with(df, points(x, fitted(m1), col = "red"))

#with(df, plot(x, p2, main = TeX(r"($P(y_k|x_k)$)"), xlab = "X (quartiles used)", ylab = ""))
#with(df, points(x, fitted(m2a), col = "red"))

with(df, plot(x, p2, main = TeX(r"($P(y_k|x_k, a_k)$)"), xlab = "X (quartiles used)", ylab = "", cex.main = 1.5))
with(df, points(x, fitted(m2b), col = "red"))

#with(df, plot(x, p2, main = TeX(r"($P(y_k|x_k)$)"), xlab = "X (deciles used)", ylab = ""))
#with(df, points(x, fitted(m3a), col = "red"))

with(df, plot(x, p2, main = TeX(r"($P(y_k|x_k, a_k)$)"), xlab = "X (deciles used)", ylab = "", cex.main = 1.5))
with(df, points(x, fitted(m3b), col = "red"))
dev.off()

# cairo_pdf(file = "figs/fig1-example2a.pdf", width = 10)
# par(mfrow=c(2,3))
# with(df, plot(x, 1/p2, main = TeX(r"($P(y_k|x_k)$)"), xlab = "X", ylab = "Y"))
# with(df, points(x, 1/fitted(m1), col = "red"))
# 
# with(df, plot(x, 1/p2, main = TeX(r"($P(y_k|x_k)$)"), xlab = "X", ylab = "Y"))
# with(df, points(x, 1/fitted(m2a), col = "red"))
# 
# with(df, plot(x, 1/p2, main = TeX(r"($P(y_k|x_k, x_k)$)"), xlab = "X", ylab = "Y"))
# with(df, points(x, 1/fitted(m2b), col = "red"))
# 
# with(df, plot(x, 1/p2, main = TeX(r"($P(y_k|x_k)$)"), xlab = "X", ylab = "Y"))
# with(df, points(x, 1/fitted(m3a), col = "red"))
# 
# with(df, plot(x, 1/p2, main = TeX(r"($P(y_k|x_k, x_k)$)"), xlab = "X", ylab = "Y"))
# with(df, points(x, 1/fitted(m3b), col = "red"))
# dev.off()
