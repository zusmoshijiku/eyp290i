library(dplyr)

#deje la proporcion min. de días con dato para aceptar el promedio de un mes
#asi solo hay 5 meses con datos faltantes: enero, febrero, marzo y abril 2001, y noviembre 2013
umbral <- 0.50

#carga de datos
mp25 <- read.csv2("mp25.csv")
names(mp25) <- c("fecha", "hora", "validados", "preliminares", "no_validados", "X")
mp25$X <- NULL


#arreglo fecha
mp25$dia <- ymd(sprintf("%06d", mp25$fecha))
range(mp25$dia)


#considera valores "validos",
#en "no validos" se eliminan los que estan bajo el minimo de los validos, que son 9 datos
mp25$valor <- coalesce(mp25$validados, mp25$preliminares, mp25$no_validados)
mp25$origen <- case_when(!is.na(mp25$validados)    ~ "validado",
                         !is.na(mp25$preliminares) ~ "preliminar",
                         !is.na(mp25$no_validados) ~ "no validado",
                         TRUE                      ~ "sin dato")
table(mp25$origen)

minimo_validado <- min(mp25$validados, na.rm = TRUE)
minimo_validado
sospechoso <- mp25$origen == "no validado" & mp25$valor < minimo_validado

summary(sospechoso) # No hay datos sospechosos

sum(is.na(mp25$valor))

#promedio mensual
mensual <- mp25 %>%
  group_by(mes = floor_date(dia, "month")) %>%
  summarise(dias_mes  = days_in_month(first(dia)),
            n_validos = sum(!is.na(valor)),
            mp25_mes  = if_else(n_validos / dias_mes >= umbral, mean(valor, na.rm = TRUE), NA_real_),
            .groups = "drop")

#Revisar que no falte ningún mes en el calendario
stopifnot(identical(mensual$mes, seq(min(mensual$mes), max(mensual$mes), by = "month")))
nrow(mensual)

#Meses que quedan sin dato
mensual %>% filter(is.na(mp25_mes))
sum(is.na(mensual$mp25_mes))

#Serie de tiempo con los meses sin dato como NA
y_na <- ts(mensual$mp25_mes, start = c(2000, 1), frequency = 12)

plot(y_na)

Y <- forecast::na.interp(y_na) # Interpola usando forecast que toma en cuenta la estacionalidad

plot(Y)

# ACF y PACF
par(mfrow = c(2, 1))
acf(Y, lag.max = 24, main = "ACF mp2.5 mensual")
pacf(Y, lag.max = 24, main = "PACF mp2.5 mensual")

par(mfrow = c(1, 1))

LSTS::periodogram(Y, plot = T)


# Probar con auto.arima
mod <- forecast::auto.arima(Y, lambda = "auto", allowdrift = T, d = 0)
mod
lmtest::coeftest(mod) # Coeficientes de la parte ARMA(2,2) no son significativos

## mod$arma guarda los ordenes en ESTE orden: p, q, P, Q, s, d, D
mod$arma
p <- mod$arma[1]
q <- mod$arma[2]
P <- mod$arma[3]
Q <- mod$arma[4]
s <- mod$arma[5]
d <- mod$arma[6]
D <- mod$arma[7]
drift <- FALSE
if(sum(names(coef(mod)) == "drift") == 1){drift = TRUE}


# Probar con Arima y quitando coeficientes no significativos
mod <- forecast::Arima(Y, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = c(rep(0, 4),rep(NA,3)))
mod
lmtest::coeftest(mod)
plot(mod)

LSTS::ts.diag(c(mod$res), lag = 24) # No se cumple blancura

# Probar agregando coeficiente ma1

mod <- forecast::Arima(Y, order = c(0,0,1), seasonal = c(2,1,0), include.drift = drift)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24) # No se cumple blancura a partir del lag 3.

# Probar agregando coeficiente ma3. Modelo: ARIMA(0,0,3)(2,1,0)[12]
q <- 3
fixed = c(NA,0,NA,NA,NA,NA)
mod <- forecast::Arima(Y, order = c(0,0,q), seasonal = c(2,1,0), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24) # Pasa la blancura hasta lag 24 (19?)


# Quizás toca agregar otros regresores






