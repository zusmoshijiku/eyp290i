library(dplyr)
library(lubridate)
library(ggplot2)
library(forecast)
library(LSTS)
library(tseries)

# Este script usa los datos diarios de la estación Parque O'Higgins para
# calcular el promedio mensual de mp2.5 desde enero 2000 hasta diciembre 2025.

# Un mes se acepta si tiene al menos el 50% de los días con datos válidos.

umbral <- 0.50

# Carga de datos
mp25 <- read.csv2("mp25.csv")
names(mp25) <- c("fecha", "hora",
                 "validados", "preliminares", "no_validados", "X")
mp25$X <- NULL


# Arreglar fecha
mp25$dia <- ymd(sprintf("%06d", mp25$fecha))
range(mp25$dia)


# considera valores "validos",
# los sospechosos son los que están en "no validados"
# y por debajo del mínimo de los validos, que son 9 datos.
mp25$valor <- coalesce(mp25$validados, mp25$preliminares, mp25$no_validados)
mp25$origen <- case_when(!is.na(mp25$validados)    ~ "validado",
                         !is.na(mp25$preliminares) ~ "preliminar",
                         !is.na(mp25$no_validados) ~ "no validado",
                         TRUE                      ~ "sin dato")
table(mp25$origen)
# no validados: 105, preliminar: 2,
# sin dato: 378, validado: 9010

minimo_validado <- min(mp25$validados, na.rm = TRUE)
minimo_validado
sospechoso <- mp25$origen == "no validado" & mp25$valor < minimo_validado

summary(sospechoso) # Hay 9 datos sospechosos
mp25$valor[sospechoso] <- NA
sum(is.na(mp25$valor)) # hay 387 días sin dato

# Promedio mensual
mensual <- mp25 %>%
  group_by(mes = floor_date(dia, "month")) %>%
  summarise(dias_mes  = days_in_month(first(dia)),
            n_validos = sum(!is.na(valor)),
            mp25_mes  = if_else(n_validos / dias_mes >= umbral,
                                mean(valor, na.rm = TRUE),
                                NA_real_), .groups = "drop")

#Revisar que no falte ningún mes en el calendario
stopifnot(identical(mensual$mes,
                    seq(min(mensual$mes), max(mensual$mes), by = "month")))
nrow(mensual)

# Meses que quedan sin dato
mensual %>% filter(is.na(mp25_mes))
sum(is.na(mensual$mp25_mes))

# Sin dato: Enero, febrero, marzo y abril 2001, y noviembre 2013.
# Todos con 0 o 1 día válido.

# Serie de tiempo con los meses sin dato como NA

y_na <- ts(mensual$mp25_mes, start = c(2000, 1), frequency = 12)
y_na <- window(y_na, start = c(2000, 1), end = c(2022, 12))

plot(y_na, main = "Serie de tiempo mensual de mp2.5 con NA",
     ylab = "mp2.5", xlab = "Año")

# largo de la serie

length(y_na) # 312 meses, serían 307 sin contar los 5 meses sin dato.

# Interpola usando forecast que toma en cuenta la estacionalidad
y <- forecast::na.interp(y_na)

plot(y, main = "Serie de tiempo mensual de mp2.5 interpolada",
     ylab = "mp2.5", xlab = "Año")

# ACF y PACF
par(mfrow = c(1, 2))
acf(y, lag.max = 24, main = "ACF mp2.5 mensual")
pacf(y, lag.max = 24, main = "PACF mp2.5 mensual")

per <- LSTS::periodogram(y, plot = TRUE)
per
max(per$lambda) # = pi

# Tendencia (esto no me funciona)
mensual$tendencia <- as.numeric(ma(y, order = 12))
ggplot(mensual, aes(x = mes)) +
  geom_line(aes(y = mp25_completa), color = "grey70", linewidth = 0.5) +
  geom_line(aes(y = tendencia), color = "darkblue", linewidth = 0.9, na.rm = TRUE) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "MP2.5 mensual y media móvil de 12 meses", x = "Fecha", y = "Concentración (µg/m³)") +
  theme_minimal()

###Ya está separado train de test
# Box-Cox
lambda <- forecast::BoxCox.lambda(y)
y.bc <- forecast::BoxCox(y, lambda = lambda)
#graficar serie con transformación de box cox lambda = 0.06
plot(y.bc, ylab = "", main = "Box-Cox(x), lambda = 0.06")

adf.test(y.bc, k = 12)
pp.test(y.bc)
kpss.test(y.bc, null = "Trend")

# efecto de la diferenciación

xs <- diff(y.bc, lag = 12)

adf.test(xs)
pp.test(xs)
kpss.test(xs, null = "Level")
kpss.test(xs, null = "Trend")
nsdiffs(xs)
ndiffs(xs, test = "adf", type = "level")

# no es necesario seguir diferenciando.



##########################################
# Probar modelos (sobre el transformado) #
##########################################

###1) AUTO ARIMA para comenzar 
mod <- forecast::auto.arima(y.bc, lambda = "auto", allowdrift = TRUE, seasonal=TRUE)
mod #ARIMA(1,0,0)(1,1,0)[12] with drift 
lmtest::coeftest(mod)
mod$arma
drift <- FALSE
if(sum(names(coef(mod)) == "drift") == 1){drift = TRUE}

##Evaluar modelo
fit <- Arima(y.bc, order = c(1,0,0), seasonal = c(1,1,0), include.drift = TRUE)
#Coeficientes significativos
lmtest::coeftest(fit)
#Raíces: estacionariedad e invertibilidad
plot(fit) #Solo tiene AR
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(fit$res), lag = 24) #hasta lag 18 hay blancura
#ACF y PACF de residuos
res <- residuals(fit)
par(mfrow = c(1, 2))
acf(res,  lag.max = 36, main = "ACF de los residuos", ylim = c(-0.3, 0.3))
pacf(res, lag.max = 36, main = "PACF de los residuos", ylim = c(-0.3, 0.3))



###2) Probar agregando coeficiente MA1
drift <- TRUE
mod <- forecast::Arima(y.bc,
  order= c(1, 0, 1), #le agrega MA 1
  seasonal=list(order = c(1, 1, 0), period = 12),
  include.drift = drift)
#Coeficientes significativos
lmtest::coeftest(mod)
#Raíces: estacionariedad e invertibilidad
plot(mod)
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(mod$residuals), lag = 24) #hasta el lag 18 hay blancura
#ACF y PACF de residuos
res <- residuals(mod)
par(mfrow = c(1, 2))
acf(res,  lag.max = 24, main = "ACF de los residuos",  ylim = c(-0.3, 0.3))
pacf(res, lag.max = 24, main = "PACF de los residuos", ylim = c(-0.3, 0.3))


###3) Probar agregando coeficiente MA seasonal: falla, las raíces se descontrolan :(
drift <- TRUE
mod <- forecast::Arima(y.bc,
                       order= c(1, 0, 0), #mantiene AR(1)
                       seasonal=list(order = c(1, 1, 1), period = 12), #agrega un MAs(1)
                       include.drift = drift)
#Coeficientes significativos
lmtest::coeftest(mod)
#Raíces: estacionariedad e invertibilidad
plot(mod)
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(mod$residuals), lag = 24) #blancura hasta 24 pero no sirve por raíces
#ACF y PACF de residuos
res <- residuals(mod)
par(mfrow = c(1, 2))
acf(res,  lag.max = 24, main = "ACF de los residuos",  ylim = c(-0.3, 0.3))
pacf(res, lag.max = 24, main = "PACF de los residuos", ylim = c(-0.3, 0.3))


###4) Probar agregando coeficiente MA1 y MA6 el resto 0 (MA2, MA3, MA4, MA5)
drift <- TRUE
mod <- forecast::Arima(y.bc,
                       order= c(1, 0, 6), #le agrega MA hasta 6
                       seasonal=list(order = c(1, 1, 0), period = 12), #mantiene esto
                       include.drift = drift,
                       fixed=c(NA, NA, 0, 0, 0, 0, NA, NA, NA)) #ar1, ma1, ma2, ma3, ma4, ma5, ma6, sar1, drift
#Coeficientes significativos
lmtest::coeftest(mod)
#Raíces: estacionariedad e invertibilidad
plot(mod)
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(mod$residuals), lag = 24) #hasta el lag 18 hay blancura
#ACF y PACF de residuos
res <- residuals(mod)
par(mfrow = c(1, 2))
acf(res,  lag.max = 24, main = "ACF de los residuos",  ylim = c(-0.3, 0.3))
pacf(res, lag.max = 24, main = "PACF de los residuos", ylim = c(-0.3, 0.3))


###5)Probar agregando coeficiente MA1 y MA19 el resto 0
drift <- TRUE
mod <- forecast::Arima(y.bc,
                       order= c(1, 0, 19), #le agrega MA hasta 19
                       seasonal=list(order = c(1, 1, 0), period = 12), #mantiene esto
                       include.drift = drift,
                       fixed=c(NA, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, NA)) #ar1, ma1, ma2, ma3, ma4, ma5,..., maX, sar1, drift
#Coeficientes significativos
lmtest::coeftest(mod)
#Raíces: estacionariedad e invertibilidad
plot(mod)
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(mod$residuals), lag = 24) #hasta el lag 23 hay blancura
#ACF y PACF de residuos
res <- residuals(mod)
par(mfrow = c(1, 2))
acf(res,  lag.max = 24, main = "ACF de los residuos",  ylim = c(-0.3, 0.3))
pacf(res, lag.max = 24, main = "PACF de los residuos", ylim = c(-0.3, 0.3))


###6) Probar agregando coeficiente AR seasonal extra
drift <- TRUE
mod <- forecast::Arima(y.bc,
                       order= c(1, 0, 0), #mantiene AR(1)
                       seasonal=list(order = c(2, 1, 0), period = 12), #agrega un MAs(1)
                       include.drift = drift,
                       fixed=c(NA, NA, NA, NA)) #ar1, sar1, sar2, drift
#Coeficientes significativos
lmtest::coeftest(mod)
#Raíces: estacionariedad e invertibilidad
plot(mod)
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(mod$residuals), lag = 24) #blancura hasta 24
#ACF y PACF de residuos
res <- residuals(mod)
par(mfrow = c(1, 2))
acf(res,  lag.max = 24, main = "ACF de los residuos",  ylim = c(-0.3, 0.3))
pacf(res, lag.max = 24, main = "PACF de los residuos", ylim = c(-0.3, 0.3))


###7) Probar agregando coeficiente AR seasonal extra y un MA19 (con intermedios 0)
drift <- TRUE
mod <- forecast::Arima(y.bc,
                       order= c(1, 0, 19), #mantiene AR(1)
                       seasonal=list(order = c(2, 1, 0), period = 12), #agrega un MAs(1)
                       include.drift = drift,
                       fixed=c(NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NA, NA, NA, NA)) #ar1, sar1, sar2, drift
#Coeficientes significativos
lmtest::coeftest(mod)
#Raíces: estacionariedad e invertibilidad
plot(mod)
#SR, ACF res y Box-Ljung
LSTS::ts.diag(c(mod$residuals), lag = 24) #blancura hasta 24
#ACF y PACF de residuos
res <- residuals(mod)
par(mfrow = c(1, 2))
acf(res,  lag.max = 24, main = "ACF de los residuos",  ylim = c(-0.3, 0.3))
pacf(res, lag.max = 24, main = "PACF de los residuos", ylim = c(-0.3, 0.3))
