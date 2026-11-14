library(dplyr)
library(lubridate)
library(ggplot2)
library(forecast)
library(MASS)
library(LSTS)
library(tseries)

#promedio mensual de mp2.5
#diario SINCA del 1 de enero de 2000 al 30 de diciembre de 2025
#el directorio tiene que ser la carpeta proyecto
#se puede setear en Session > Set Working Directory > To Source File Location

#########################

# carga y limpieza de datos

#########################

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
mp25 %>% filter(sospechoso) %>% select(dia, valor, origen)
mp25$valor[sospechoso] <- NA
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
x_na <- ts(mensual$mp25_mes, start = c(2000, 1), frequency = 12)

#se rellenar los meses sin dato con interpolación estacional
#interpola linealmente la serie sin ese ciclo y después se lo vuelve a sumar
x <- na.interp(x_na)
mensual$estimado <- is.na(mensual$mp25_mes)
mensual$mp25_completa <- as.numeric(x)

#valores estimados(los 9 que se borraron)
mensual %>% filter(estimado) %>% select(mes, n_validos, mp25_completa)

#serie lista
x

#########################

# analisis exploratorio

#########################

#gráfico de la serie en el tiempo
ggplot(mensual, aes(x = mes, y = mp25_completa)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title="MP 2.5 promedio mensual",x = "Fecha", y = "Concentración (µg/m³)") +
  theme_minimal()

# ACF y PACF
par(mfrow = c(2, 1))
acf(x, lag.max = 24, main = "ACF mp2.5 mensual")
pacf(x, lag.max = 24, main = "PACF mp2.5 mensual")

par(mfrow = c(1, 1))


#periodograma
per <- LSTS::periodogram(x, plot = FALSE)
periodograma <- data.frame(lambda = per$lambda, I = per$periodogram) %>%
  mutate(periodo_meses = 2 * pi / lambda)
ggplot(periodograma, aes(x = lambda, y = I)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  labs(title = "Periodograma mp2.5 mensual", x = "frecuencia", y = "") +
  theme_minimal()

periodograma %>% arrange(desc(I)) %>% head(5)

#tendencia
mensual$tendencia <- as.numeric(ma(x, order = 12))
ggplot(mensual, aes(x = mes)) +
  geom_line(aes(y = mp25_completa), color = "grey70", linewidth = 0.5) +
  geom_line(aes(y = tendencia), color = "darkblue", linewidth = 0.9, na.rm = TRUE) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "MP2.5 mensual y media móvil de 12 meses", x = "Fecha", y = "Concentración (µg/m³)") +
  theme_minimal()

#estacionalidad
ggplot(mensual, aes(x = factor(month(mes)), y = mp25_completa)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "MP2.5 por mes del año", x = "Mes", y = "Concentración (µg/m³)") +
  theme_minimal()

#atípicos
atipicos <- tsoutliers(x)
data.frame(mes = mensual$mes[atipicos$index],
           observado = as.numeric(x)[atipicos$index],
           esperado  = atipicos$replacements)


###########################
# Probar modelo con Box-Cox
###########################

##Box Cox lambda = 0.24
lambda <- forecast::BoxCox.lambda(x) 
x.bc <- forecast::BoxCox(x, lambda = lambda)
#graficar serie con transformación de box cox
plot(x.bc, ylab = "", main = "Box-Cox(x), lambda = 0.24")

#ADF, PP y KPSS para ver estacionariedad
adf.test(x.bc)
pp.test(x.bc) 
kpss.test(x.bc) 

#ver diferenciación
D <- forecast::nsdiffs(x.bc) #D=1
d <- forecast::ndiffs(diff(x.bc, lag = 12)) #d=0

#revisar la diferenciación: todos de acuerdo con que no más
xs <- diff(x.bc, lag = 12)
ndiffs(xs, test = "kpss")
ndiffs(xs, test = "adf")
ndiffs(xs, test = "pp")

#ACF y PACF sobre la serie transformada con Box Cox y 1 dif seasonal
par(mfrow = c(2, 1))
acf(xs, lag.max = 24, main = "ACF BC diff seasonal")
pacf(xs, lag.max = 24, main = "PACF BC diff seasonal")

#Ver que modelo sugiere auto arima seasonal
fit_auto <- auto.arima(x.bc, seasonal = TRUE) #aqui se usa la sin diferenciar pues autoarima decide
summary(fit_auto) #sugiere un ARIMA(1,0,1)(0,1,1)[12]

#Probar el modelo
fit <- Arima(x.bc, order = c(1,0,1), seasonal = c(0,1,1), include.drift = TRUE)
checkresiduals(fit) #Box-Ljung no rechaza la blancura

#ver los residuos
res <- residuals(fit)

#ACF y PACF de los residuos
par(mfrow = c(2, 1))
acf(res,  lag.max = 36, main = "ACF de los residuos (no B-C)")
pacf(res, lag.max = 36, main = "PACF de los residuos (no B-C)")
par(mfrow = c(1, 1))


###########################
# Probar modelo sin Box-Cox
###########################

#ADF, PP y KPSS para ver estacionariedad
adf.test(x)
pp.test(x) 
kpss.test(x) 

#ver diferenciación
D <- forecast::nsdiffs(x) #D=1
d <- forecast::ndiffs(diff(x, lag = 12)) #d=0

#revisar la diferenciación: todos de acuerdo con que no más
xs <- diff(x, lag = 12)
ndiffs(xs, test = "kpss")
ndiffs(xs, test = "adf")
ndiffs(xs, test = "pp")

#ACF y PACF sobre la serie 1 dif seasonal
par(mfrow = c(2, 1))
acf(xs, lag.max = 24, main = "ACF BC diff seasonal")
pacf(xs, lag.max = 24, main = "PACF BC diff seasonal")

#Ver que modelo sugiere auto arima seasonal
fit_auto <- auto.arima(x, seasonal = TRUE) #aqui se usa la sin diferenciar pues autoarima decide
summary(fit_auto) #sugiere un ARIMA(1,0,1)(1,1,2)[12]

#Probar el modelo
fit <- Arima(x, order = c(1,0,1), seasonal = c(1,1,2), include.drift = TRUE)
checkresiduals(fit) #Box-Ljung rechaza la blancura :(

#ver los residuos
res <- residuals(fit)

par(mfrow = c(2, 1))
acf(res,  lag.max = 36, main = "ACF de los residuos (no B-C)")
pacf(res, lag.max = 36, main = "PACF de los residuos (no B-C)")
par(mfrow = c(1, 1))