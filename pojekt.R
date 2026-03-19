library(quantmod)
nazwy_cen <- c("BTC-USD", "ETH-USD", "XRP-USD","^DJI")
prices <- xts()
for (stock_index in 1:length(nazwy_cen))
  prices <- cbind(prices, Ad(getSymbols(nazwy_cen[stock_index], 
                                        from = "2021-01-01", to = "2024-12-31", auto.assign = FALSE)))
colnames(prices) <- nazwy_cen
indexClass(prices) <- "Date"
#View(prices)
#head(prices)

# plot the three series of log-prices
plot(log(prices), col = c("brown", "blue", "magenta","green"),
     main = "Logarytmy cen", legend.loc = "left")