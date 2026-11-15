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

plot(y_na, main = "Serie de tiempo mensual de mp2.5 con NA",
     ylab = "mp2.5", xlab = "Año")

# largo de la serie

length(y_na) # 312 meses, serían 307 sin contar los 5 meses sin dato.

# Interpola usando forecast que toma en cuenta la estacionalidad
y <- forecast::na.interp(y_na)

plot(y, main = "Serie de tiempo mensual de mp2.5 interpolada",
     ylab = "mp2.5", xlab = "Año")

# ACF y PACF
par(mfrow = c(2, 1))
acf(y, lag.max = 24, main = "ACF mp2.5 mensual")
pacf(y, lag.max = 24, main = "PACF mp2.5 mensual")

par(mfrow = c(1, 1))

per <- LSTS::periodogram(y, plot = TRUE)

max(per$lambda) # = pi

# Tendencia
mensual$tendencia <- as.numeric(ma(x, order = 12))
ggplot(mensual, aes(x = mes)) +
  geom_line(aes(y = mp25_completa), color = "grey70", linewidth = 0.5) +
  geom_line(aes(y = tendencia), color = "darkblue", linewidth = 0.9, na.rm = TRUE) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "MP2.5 mensual y media móvil de 12 meses", x = "Fecha", y = "Concentración (µg/m³)") +
  theme_minimal()

# Box-Cox

lambda <- forecast::BoxCox.lambda(x)
y.bc <- forecast::BoxCox(x, lambda = lambda)
#graficar serie con transformación de box cox lambda = 0.24
plot(y.bc, ylab = "", main = "Box-Cox(x), lambda = 0.24")

adf.test(y.bc)
pp.test(y.bc)
kpss.test(y.bc)
