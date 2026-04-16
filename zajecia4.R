library(rmgarch)
data(dji30retw)

#View(dji30retw)
Dat <-dji30retw[, 1:3]
y<-Dat#na tym wektorze pracujemy

dim(Dat)
names(y)

specyfikacja <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(1,0), include.mean = FALSE), distribution.model="norm")

dcc_specyfikacja<-dccspec(uspec = multispec( replicate(3, specyfikacja)), VAR = TRUE, lag =1, dccOrder = c(1,1), distribution = "mvnorm")

fit_dcc <- dccfit(spec = dcc_specyfikacja, data = y)

?dccfit

names(fit_dcc@mfit)

class(fit_dcc)

slotNames(fit_dcc)
names(fit_dcc@model)
names(fit_dcc@mfit)

fit_dcc@model$varcoef

names(y)

plot(fit_dcc)

simulation<-dccsim(fit_dcc,n.sim=3000,n.start=1000,m.sim=1, rseed=1)
slotNames(simulation)
names(simulation@msim)
(simulation@msim$simQ)
simulation@msim$rseed
names(simulation@msim)
(simulation@msim$simX)
simulation@msim$simH[[1]][,,111] #Generalnei to jest warunkowa kowaranacja w czasie - 111 - oznacza czas
simulation@msim$simH[[1]][1,1,111]
simulation@msim$simR[[1]][,,111] #Macierz warunkowych korelacji
simulation@msim$simR[[1]][1,2,111]


#Wykres wygenerowanych szeregow z modelu DDC
par(mfrow=c(3,1))
plot(simulation@msim$simX[[1]][,1],type="l")#simulated time series with DCC-GARCH process
plot(simulation@msim$simX[[1]][,2],type="l")
plot(simulation@msim$simX[[1]][,3],type="l")

# Wykresy warunkowych wariancji dla wygenerowanych szeregow
par(mfrow=c(3,1))
plot(simulation@msim$simH[[1]][1,1,],type="l")#simulated time series with DCC-GARCH process
plot(simulation@msim$simH[[1]][2,2,],type="l")
plot(simulation@msim$simH[[1]][3,3,],type="l")


# Wykresy korelacji dla wygenerowanych szeregow - kazdy wykres z kazdym trzeba skorelowac
par(mfrow=c(3,1))
plot(simulation@msim$simR[[1]][1,2,],type="l")
plot(simulation@msim$simR[[1]][2,3,],type="l")
plot(simulation@msim$simR[[1]][3,1,],type="l")


cluster<-NULL

Dat = dji30retw[, 1:3, drop = FALSE]

uspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                   distribution.model = "norm")
spec1 = cgarchspec(uspec = multispec( replicate(3, uspec) ), VAR = TRUE, robust = FALSE, lag = 0, lag.max = NULL,
                   lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL,
                   robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500),
                   dccOrder = c(1,1), asymmetric = FALSE, distribution.model = list(copula = c("mvnorm", "mvt")[1],
                                                                                    method = c("Kendall", "ML")[2], time.varying = FALSE,
                                                                                    transformation = c("parametric", "empirical", "spd")[1]))
#CCC-MGARCH
fit1 = cgarchfit(spec1, data = Dat, cluster = cluster, fit.control = list(eval.se=FALSE))
(fit1)
#Teraz generujemy z modelu CCC-MGARCH, który wcześniej został dopasowany do rzeczywistych obserwacji, dzięki temu mamy realne wartości parametrów
simulationCCC = cgarchsim(fit1, n.sim = 1000, m.sim = 1, startMethod = "sample", cluster=cluster)

#Do danych z CCC dopasuje DCC

specyfikacja = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model = list(armaOrder = c(1,0), include.mean = FALSE), distribution.model="norm")
dcc_specyfikacja<-dccspec(uspec = multispec( replicate(3, specyfikacja)), VAR = TRUE, lag =0, dccOrder = c(1,1), distribution = "mvnorm")
fit_dcc_symulacja<-dccfit(spec = dcc_specyfikacja, data = simulationCCC@msim$simX[[1]])

#generujemy dane z modelu DCC:
simulation<-dccsim(fit_dcc_symulacja,n.sim=3000,n.start=1000,m.sim=1, rseed=1)



#Wykresy wygenerowanych szeregów z modelu DCC
par(mfrow=c(3,1))
plot(simulation@msim$simX[[1]][,1],type="l")#simulated time series with DCC-GARCH process
plot(simulation@msim$simX[[1]][,2],type="l")
plot(simulation@msim$simX[[1]][,3],type="l")


#Wykresy wygenerowanych szeregów z modelu DCC z ograniczeniem na osi Y 
par(mfrow=c(3,1))
plot(simulation@msim$simR[[1]][1,2,],type="l",ylim=c(0.2,0.6))
abline(h=median(simulation@msim$simR[[1]][1,2,]),col='red')
plot(simulation@msim$simR[[1]][1,3,],type="l",ylim=c(0.2,0.6))
abline(h=median(simulation@msim$simR[[1]][1,3,]),col='red')
plot(simulation@msim$simR[[1]][2,3,],type="l",ylim=c(0.2,0.6))
abline(h=median(simulation@msim$simR[[1]][2,3,]),col='red')



simulationCCC <- cgarchsim(fit1, n.sim = 3000, m.sim = 1, startMethod = "sample", cluster=cluster)

y<-simulationCCC@msim$simX[[1]][,1:3]

DCCtest(Dat, garchOrder = c(1,1), n.lags = 2, solver = "solnp", solver.control = list(), cluster = NULL, Z = NULL)

simulationDCC<-dccsim(?)

y<-simulationDCC@msim?
  
DCCtest(?, garchOrder = ?, n.lags = ?, solver = "solnp", solver.control = list(), cluster = NULL, Z = NULL)

spec<-ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=TRUE, garchInMean = FALSE, inMeanType = 2), distribution.model="sstd", fixed.pars=list(mu=0.001,omega=0.00001, alpha1=0.05, beta1=0.90, shape=4,skew=2))

y<-ugarchpath(spec, n.sim=3000, n.start=1000, m.sim=3)



fit1@mfit$convergence


simulationCCC <- cgarchsim(fit1, n.sim = 3000, m.sim = 1, startMethod = "sample", cluster=cluster,rseed=1)
y<-simulationCCC@msim$simX[[1]][,1:3]
DCCtest(y, garchOrder = c(1,1), n.lags = 2, solver = "solnp",
        solver.control = list(), cluster = NULL, Z = NULL)


simulationDCC<-dccsim( fit_dcc_symulacja@mfit$coef,n.sim=3000,n.start=1000,m.sim=1, rseed=1)
y<-simulationDCC@msim$simX[[1]][,1:3]
DCCtest(y, garchOrder = c(1,1), n.lags = 2, solver = "solnp",
        solver.control = list(), cluster = NULL, Z = NULL)

simulationDCC<-dccsim(fit_dcc_symulacja,n.sim=3000,n.start=1000,m.sim=1, rseed=1)

fit_dcc_symulacja@mfit$coef[13]=0.4
fit_dcc_symulacja@mfit$coef[14]=0.2
fit_dcc_symulacja@mfit$coef
fit_dcc_symulacja@mfit$matcoef
fit_dcc_symulacja@mfit$matcoef[13]=0.4
fit_dcc_symulacja@mfit$matcoef[14]=0.2
fit_dcc_symulacja

simulationDCC<-dccsim(fit_dcc_symulacja,n.sim=3000,n.start=1000,m.sim=1, rseed=1)
y<-simulationDCC@msim$simX[[1]][,1:3]
DCCtest(y, garchOrder = c(1,1), n.lags = 2, solver = "solnp",
        solver.control = list(), cluster = NULL, Z = NULL)

simulationDCC<-dccsim(fit_dcc,n.sim=3000,n.start=1000,m.sim=1, rseed=9)
y<-simulationDCC@msim$simX[[1]][,1:3]
DCCtest(y, garchOrder = c(1,1), n.lags = 2, solver = "solnp",
        solver.control = list(), cluster = NULL, Z = NULL)

plot(simulation@msim$simR[[1]][1,2,],type="l",ylim=c(0.2,0.6))
plot(y)
plot(DCCtest)
# WYKRES DLA WARUNKOWYCH KORELACJI
plot(simulation@msim$simR[[1]][1,2,],type="l",ylim=c(-1,1))

# WYKRES DLA WARUNKOWYCH KORELACJI PRZY ODRZUCENIU HIPOTEZY 0

#########
plot(fit_dcc, which = 4, series = c(1, 3))
plot(fit_dcc, which = 4, series = c(2, 3))


# Zadanie 5 do domu
