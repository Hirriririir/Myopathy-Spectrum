library(DescTools)
library(readxl)

# https://andrisignorell.github.io/DescTools/reference/JonckheereTerpstraTest.html

g <- ordered(rep(1:5, rep(10,5)))
x <- rnorm(50) + 0.3 * as.numeric(g)

JonckheereTerpstraTest(x, g)

# CDM #----
CDM <- read_excel("Validation/Jonckheere trend test.xlsx", sheet = "CDM")
CDM$X_category <- ordered(CDM$X_category)

JonckheereTerpstraTest(CDM$CTG_Expansion, CDM$X_category)

# LGMD R12 #----
LGMDR12 <- read_excel("Validation/Jonckheere trend test.xlsx", sheet = "LGMDR12")
LGMDR12$X_category <- ordered(LGMDR12$X_category)

JonckheereTerpstraTest(LGMDR12$Mercuri_score, LGMDR12$X_category)

JonckheereTerpstraTest(LGMDR12$`6MWD`, LGMDR12$X_category)

JonckheereTerpstraTest(LGMDR12$`10MWT`, LGMDR12$X_category)


# FSHD #----
FSHD <- read_excel("Validation/Jonckheere trend test.xlsx", sheet = "FSHD")
FSHD$X_category <- ordered(FSHD$X_category)

JonckheereTerpstraTest(FSHD$Path_score, FSHD$X_category)

JonckheereTerpstraTest(FSHD$CSS, FSHD$X_category)

JonckheereTerpstraTest(FSHD$Fat_fraction, FSHD$X_category)