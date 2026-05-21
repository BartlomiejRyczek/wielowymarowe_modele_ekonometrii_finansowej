dane1<- read.csv2(file="pkn_20152021.csv",dec=".",sep=",")

dane1<-as.data.frame(dane1[,c(1,5)])
names(dane1)[2]<-"pkn"
dane2<-read.csv(file="kgh_20152021.csv", sep=",", dec=".")
dane2<-as.data.frame(dane2[,c(1,5)])
names(dane2)[2]<-"kgh"

dane3<-read.csv(file="pko_d.csv", dec=".")
dane3<-as.data.frame(dane3[,c(1,5)])
names(dane3)[2]<-"PKOBP"


lnrdane1 <- diff(log(dane1$pkn))
lnrdane2 <- diff(log(dane2$kgh))
lnrdane3 <- diff(log(dane3$PKOBP))

y<-na.omit(merge(dane1,dane2, by="Data"))
y<-na.omit(merge(y,dane3, by="Data"))

y<-y[,c(2,3,4)]
y_returns<-y[-1,]
y_returns[,1]<-diff(log(y[,1]))
y_returns[,2]<-diff(log(y[,2]))
y_returns[,3]<-diff(log(y[,3]))
nobs<-dim(y_returns)[1]
plot(y_returns[,3],type="l")

# nie robimy dopasowania diagnostyki reszt bo nie ma czasu 
install.packages("MTS")
install.packages("rmgarch")
install.packages("vars")
library(vars)
library(rmgarch)
library(MTS)

varx_fit_returns<-varxfit(y_returns,p=2, postpad = "constant")
yspec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(2,2)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution.model="sstd")#Tu False bo b?d? mia? VAR zamist constance
yspec = multispec(replicate(3, yspec1))
spec_dcc<-dccspec(uspec = yspec, VAR = TRUE, lag =2, dccOrder = c(1,1), distribution = "mvt")
multf = multifit(yspec, data = varx_fit_returns$xresiduals,cluster=NULL)
yfit_DCC_returns<-dccfit(spec_dcc, data = y_returns, fit.control = list(eval.se=TRUE),VAR.fit = varx_fit_returns, fit = multf, cluster = NULL)

MarchTest(yfit_DCC_returns@mfit$stdresid,lag=4)#jeden z testów odrzucil tj. brak efekt ARCH w resztach standaryzowanych
arch.test(VAR(yfit_DCC_returns@mfit$stdresid))#standaryzowane reszty

arch.test(VAR(yfit_DCC_returns@model$residuals))#zwykle reszty w sumie nie potrzebny bo testujemy na standaryzowanych resztach

K<-3# 3 walory ryzykowne
H_t<-matrix(0,nrow=K,ncol=K)#macierz war. kowariancji dla ust t
C_t<-matrix(0,nrow=1,ncol=1)#stala normujaca dla chwili t
w_wagi<-matrix(0,nrow=(nobs+1),ncol=K)
w_wagi_statyczny<-matrix(0,nrow=(nobs+1),ncol=K)+c(1/3,1/3,1/3)
i_jedn<-1+vector("integer",length=K)
alfa<-0.05#poziom istotności
VaR_Long_returns<-matrix(0,nrow=(nobs))
VaR_Long_portfolio<-matrix(0,nrow=(nobs))
est_mi<-matrix(0,nrow=nobs,ncol=K)#Prognozy warunkowych wartosci oczekiwanych; Uzyskujemy z dopasowania modelu VAR
est_mi<-varx_fit_returns$xfitted#wartosci predyktywne = srednie =oszacowania modelu

V_assets<-matrix(0,nrow=(nobs+1),ncol=K)#+1 bo uwzgledniam warunek startowy; Wartosci instrument tworzacych portfel
V_0_assets<-as.matrix(y[1,c(1,2,3)])#Poczatkowe wartosci instrumentow; 
V_assets[1,]<-V_0_assets#Poczatkowe wartosci instrumentow; tutaj arbitralnie przyjete
V_assets<-as.matrix(y[,c(1,2,3)])


par(mfrow = c(3, 1)) 

plot(V_assets[, 1], type = "l", main = "PKN", ylab = "Wartość", xlab = "Czas")
plot(V_assets[, 2], type = "l", main = "KGH", ylab = "Wartość", xlab = "Czas")
plot(V_assets[, 3], type = "l", main = "PKOBP", ylab = "Wartość", xlab = "Czas")

par(mfrow = c(1, 1))

V_portfolio<-vector("integer",length=(nobs+1))#Values of the portfolio
V_portfolio[1]<-1000#arbitralnie przyjeta warto?? portfela
V_portfolio_statyczny<-vector("integer",length=(nobs+1))#Values of the portfolio
V_portfolio_statyczny[1]<-1000#arbitralnie przyjeta wartosc portfela


# W statycznym portfelu inwestujemy w akcje po tyle samo czyli po 1/3 w naszym przypadku


w_wagi[1,]<-c(1/3,1/3,1/3)
liczbaAkcji<-matrix(0,nrow=(nobs+1),ncol=K)
liczbaAkcji[1,]<-V_portfolio[1]*w_wagi[1,]/V_assets[1,]#liczba akcji odpowiadajaca momentowi startowemu 1
liczbaAkcji_statyczny<-matrix(0,nrow=(nobs+1),ncol=K)
liczbaAkcji_statyczny[1,]<-V_portfolio_statyczny[1]*w_wagi_statyczny[1,]/V_assets[1,]

alfa<-0.05

for(i in 2:(nobs+1)){
  H_t<-yfit_DCC_returns@mfit$H[,,i-1]
  C_t<-(t(as.matrix(i_jedn))%*%solve(H_t))%*%as.matrix(i_jedn)
  w_wagi[i,]<-solve(H_t)%*%as.matrix(i_jedn/C_t[1,1]) 
  VaR_Long_returns[i-1]<--sum(w_wagi[i-1,]*est_mi[i-1,])-qnorm(alfa)*((t(as.matrix(w_wagi[i-1,]))%*%H_t)%*%as.matrix(w_wagi[i-1,]))^0.5#Wz?r poprawny przy za?o?eniu warunkowego rozk?adu normalnego ale jak si? okazuje dzia?a
  V_portfolio[i]<-sum(liczbaAkcji[i-1,]*(V_assets[i,]))#sum(w_wagi[i,]*(V_assets[i,]))
  V_portfolio_statyczny[i]<-sum(liczbaAkcji_statyczny[i-1,]*(V_assets[i,]))
  liczbaAkcji[i,]<-V_portfolio[i]*w_wagi[i,]/V_assets[i,]
  liczbaAkcji_statyczny[i,]<-V_portfolio_statyczny[i]*w_wagi_statyczny[i,]/V_assets[i,]
  VaR_Long_portfolio[i-1]<-V_portfolio[i-1]*(exp(sum(w_wagi[i-1,]*est_mi[i-1,])+qnorm(alfa)*((t(as.matrix(w_wagi[i-1,]))%*%H_t)%*%as.matrix(w_wagi[i-1,]))^0.5))
}


par(mfrow=c(3,2))
plot(w_wagi[, 1], type = "l")
plot(w_wagi[, 2], type = "l")
plot(w_wagi[, 3], type = "l")

plot(w_wagi_statyczny[, 1], type = "l")
plot(w_wagi_statyczny[, 2], type = "l")
plot(w_wagi_statyczny[, 3], type = "l")



portfolio_returns<-diff(log(V_portfolio))
portfolio_returns_statyczny<-diff(log(V_portfolio_statyczny))


par(mfrow=c(1,1))
#pozycja dluga to na dole pozycja krotka to gorna linia 
matplot(portfolio_returns,type="l")
matplot(-VaR_Long_returns,col="red",add=T,type="l")
matplot(VaR_Long_returns,col="red",add=T,type="l")

#narysuj wartosc portflea statycznego czerwony i dynamicznego zielony na jednym wykresie
#kod z chatu
par(mfrow = c(1, 1))

plot(V_portfolio_statyczny, type = "l",col = "red",lwd = 2,
     main = "Wartość portfela statycznego i dynamicznego",xlab = "Czas",ylab = "Wartość portfela",
     ylim = range(c(V_portfolio_statyczny, V_portfolio), na.rm = TRUE))

lines(V_portfolio, col = "green", lwd = 2)
#kod z zajec 
par(mfrow=c(1,1))
matplot(V_portfolio, col = "darkgreen", type = "l")
matplot(V_portfolio_statyczny,col="red",add=T,type="l")
#wnioski skala ma znaczenie

library("rugarch")
print(VaRTest(alfa, portfolio_returns, -VaR_Long_returns))
print(VaRTest(alfa, V_portfolio[-1], VaR_Long_portfolio))
#interpretacja
# nie ma podstaw do odrzcenia hipoezy 0 ze przekroczenia wyosza 5%


