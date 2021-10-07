library(DescTools)
library(wavScalogram)

# PHENOLOGICAL DIVERSITY

# Description

# qPD is a measure of phenological diversity based on the framework proposed by Chao, et al. (2014).

# This tool calculates a phenological diversity measure from the Hill numbers perspective. First, the script transforms data on discrete phenological species patterns into modelled continuous patterns through wavelet analysis. The user can modify the "tau" parameter and time distance of sampling effort.
# Second, phenological diversity is calculated for these modelled phenological patterns, where the "q" parameter controls the weight species phenological intensity has on the diversity measure.

## Data preparation throught wavelets transform ##

# Time distance between samples of studied phenological trait
dt <- 1 # 1 year between each phenological measurement 

# Parameter that controls the wiggliness of the phenological pattern: t = 2, default, t -> 0 smoothed slope, t -> \infty  extremely wiggly. Depending on the data set at hand, the user may change this default: wiggly when the phenological pattern changes abruptly over time.
tau <- 2

# SECTION 1: Data transformation from discrete species phenological patterns into continuous modelled ones.

# Intensity data must be arranged as: colunms -> species; rows -> sampling time.
setwd("...")
signals <- read.csv("betadata1.csv")

i = 1

n.time <- nrow(signals)
n.spp <- ncol(signals)
time <- 1:n.time
RelSpPhCumInt <- rep(NA, n.spp) # Relative Species Phenological Cumulative Intensity vector
wSpPhInt <- SpPhInt <- matrix(NA, nrow = n.time, ncol = n.spp) # Weighted and absolute Species Phenological Intensity Continuous Patterns
for (i in 1:n.spp) {
  SpPhPat <- signals[, i]
  RelSpPhCumInt[i] <- sum(SpPhPat)
  if (sd(SpPhPat) == 0)
    SpPhInt[, i] <- SpPhPat/sum(SpPhPat)
  if (sd(SpPhPat) > 0) {
    wavelet <- cwt_wst(signal = SpPhPat, dt = dt, makefigure = FALSE) # wavelet analysis
    z <- abs(wavelet$coefs)^tau # wavelet-transformed signal
    SpPhInt[, i] <- rowSums(z)/sum(z) # modelled phenological continuous pattern
  }
}
RelSpPhCumInt <- RelSpPhCumInt/sum(RelSpPhCumInt)
for(i in 1:n.spp)
  wSpPhInt[, i] <- RelSpPhCumInt[i]*SpPhInt[, i]

# Plotting the community phenological patterns
for(i in 1:ncol(signals)) {
  if (i == 1)
    plot (time, wSpPhInt[, i], type = "l", las = 1, bty = "l", xlab = "Time (yr)", ylab = "Intensity", ylim = c(0, max(wSpPhInt)))
  if (i != 1)
    lines(time, wSpPhInt[,i])
}
  
# SECTION 2: Phenological diversity calculation

q <- c(0, 1, 2) # these q-values can be changed if needed
qPD <- rep(NA, ncol = length(q)) 

Qij <- Oij <- matrix(NA, ncol = n.spp, nrow = n.spp)
for (i in 1:n.spp) {
  for (j in 1:n.spp) {
    Ii <- SpPhInt[, i]
    Ij <- SpPhInt[, j]
    Oij[i, j] <- 1-2*sum(Ii*Ij)/(sum(Ii^2) + sum(Ij^2))
    Qij[i, j] <- Oij[i, j]*RelSpPhCumInt[i]*RelSpPhCumInt[j]
  }
}
Q <- sum(Qij)
qPDij <- array(NA, dim = c(length(q), n.spp, n.spp))
for (k in 1:length(q)) {
  if (q[k] != 1 & Q != 0) {
    for (i in 1:n.spp)
      for (j in 1:n.spp)
        qPDij[k, i, j] <- Oij[i, j]*(RelSpPhCumInt[i]*RelSpPhCumInt[j]/Q)^q[k]
    qPD[k] <- (sum(qPDij[k, ,])^(1/(1-q[k])))^(0.5)
  }
  if (q[k] == 0 & Q == 0)
    qPD[k] <- 1
  if (q[k] == 2 & Q == 0) {
    for (i in 1:n.spp)
      for (j in 1:n.spp) {
        qPDij[k, i, j] <- (RelSpPhCumInt[i]*RelSpPhCumInt[j])^q[k]
      } 
    qPD[k] <- ((sum(qPDij[k, ,])^(1/(1-q[k])))^(0.5))/n.spp
  }
  if (q[k] == 1) {
    for (i in 1:n.spp) 
      for (j in 1:n.spp)
        qPDij[k, i, j] <- Oij[i, j]*RelSpPhCumInt[i]*RelSpPhCumInt[j]/Q*log(RelSpPhCumInt[i]*RelSpPhCumInt[j]/Q)
      qPD[k] <- (exp(-sum(qPDij[k, ,])))^(0.5)
      if (Q == 0) {
        for (i in 1:n.spp) 
          for (j in 1:n.spp)
            qPDij[k, i, j] <- RelSpPhCumInt[i]*RelSpPhCumInt[j]*log(RelSpPhCumInt[i]*RelSpPhCumInt[j])
        qPD[k] <- ((exp(-sum(qPDij[k, ,])))^(0.5))/n.spp
      }
    }
}

# References
# Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species qPDersity, Phylogenetic qPDersity, Functional qPDersity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297-324. <doi:10.1146/annurevecolsys-120213-091540>.