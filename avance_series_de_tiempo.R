library(dplyr)
library(ggplot2)
library(lubridate)
library(zoo)
library(forecast)
library(lmtest)

#Descargué el archivo de datos diarios de MP 2.5 desde 1 de enero de 2000 hasta 31 de diciembre de 2025 desde SINCA 
#Lo tengo guardado como mp25.csv, para cargarlo uso import dataset, from text (base)

#Nombres a columnas y elimina la última
names(mp25) <- c("fecha", "hora", "validados", "preliminares", "no_validados", "X")
mp25$X <- NULL

#Convertir fecha (rellenar con ceros a 6 dígitos) y hora (a 4 dígitos), luego cambia a col fecha_hora con formato date
mp25$fecha <- sprintf("%06d", mp25$fecha)
mp25$hora  <- sprintf("%04d", mp25$hora)
mp25$fecha_hora <- ymd(mp25$fecha) + hours(as.integer(substr(mp25$hora, 1, 2)))

#Unir validados, preeliminares y no validados en 1 col nueva llamada valor
mp25$valor <- coalesce(mp25$validados, mp25$preliminares, mp25$no_validados)

#Rangos de fecha 
range(mp25$fecha) #2000 enero 1 a 2025 diciembre 31
summary(mp25$valor) #hay 378 NA

#graficar datos diarios para ver
ggplot(mp25, aes(x = fecha_hora, y = valor)) +
  geom_line(color = "darkblue", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "MP2.5", x = "Fecha", y = "Concentración (µg/m³)") +
  theme_minimal()

#Revisar largo de intervalos faltantes (esta me la dio la IA)
intervalos_faltantes <- function(df, col_tiempo = "fecha_hora", col_valor = "valor") {
  t <- df[[col_tiempo]]
  v <- df[[col_valor]]
  es_na <- is.na(v)
  r <- rle(es_na)
  fin    <- cumsum(r$lengths)
  inicio <- fin - r$lengths + 1
  
  tibble(
    inicio_idx = inicio,
    fin_idx    = fin,
    es_na      = r$values,
    largo      = r$lengths
  ) %>%
    filter(es_na) %>%
    transmute(
      inicio = t[inicio_idx],
      fin    = t[fin_idx],
      largo_dias = largo
    )
}

faltantes <- intervalos_faltantes(mp25)
faltantes
summary(faltantes$largo_dias)
table(faltantes$largo_dias)
#En resumen:
#51 casos de 1 día faltante
#23 casos de 2 días faltantes
#4 casos de 3 días faltantes
#podemos hacer interpolación lineal para todos estos

#18 casos de entre 4 y 11 días faltantes: son bastantes, quedan como NA

#1 caso de 45 días faltantes (entre 2013-octubre-28 y 2013-diciembre-11): dejar como semanas NA
#1 caso de 119 días faltantes (2001-enero-1 y 2001-abril-29): eliminar pues es al principio
faltantes %>% filter(largo_dias %in% c(45, 119)) #este es para ver las fechas de los vacíos grandes de antes

#Filtrar para comenzar despúes del vacío largo
mp25$dia <- ymd(mp25$fecha)
mp25d <- mp25[order(mp25$dia), c("dia", "valor")]
mp25d <- mp25d[mp25d$dia >= "2001-04-30", ] #comienza después del vacío largo

#Ver los nuevos faltantes desde 30 de abril de 2001 en adelante
faltantes_d <- intervalos_faltantes(mp25d)
faltantes_d
summary(faltantes_d$largo_dias)
table(faltantes_d$largo_dias)


#Interpolar solo vacíos de 1, 2 o 3 días
mp25d$valor_interpolado <- na.approx(mp25d$valor, maxgap = 3, na.rm = FALSE)

#Serie de promedio semanal (pero solo si hay más de 4 días)
mp25d$semana <- floor_date(mp25d$dia, "week", week_start = 1)
semanal <- mp25d %>%
  group_by(semana) %>%
  summarise(n_validos = sum(!is.na(valor_interpolado)),
    mp25_sem  = if_else(n_validos >= 4, mean(valor_interpolado, na.rm = TRUE), NA_real_),
    .groups = "drop")

#graficar semanal para ver (tiene 1288 filas)
ggplot(semanal, aes(x = semana, y = mp25_sem)) +
  geom_line(color = "darkgreen", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title = "MP2.5 semanal", x = "Fecha", y = "Concentración (µg/m³)") +
  theme_minimal()

#faltantes de la serie semanal
sum(is.na(semanal$mp25_sem)) #hay 22 faltantes ahora





#Análisis de la serie

#crear la serie de tiempo
x <- ts(semanal$mp25_sem, frequency = 52)

#diferenciar?
D <- forecast::nsdiffs(x) #requiere seasonal
forecast::ndiffs(diff(x, lag = 52)) #no requiere regular después de seasonal
x_dif<-diff(x) #serie diferenciada

#graficar serie diferenciada
autoplot(x_dif) +
  labs(title = "MP2.5 semanal (diferencia regular)",
       y = "Concentración MP 2.5 (µg/m³)", x = "Año") +
  theme_minimal()

#ACF y PACF (falla pues hay NA)
par(mfrow = c(2, 1))
acf(x_dif, lag.max = 40, main = "ACF")
pacf(x_dif, lag.max = 40, main = "PACF")

#Probar modelo (solo para ver si es posible con los NA, pero faltan más cosas antes)
fit_auto <- auto.arima(x, seasonal = TRUE)
summary(fit_auto)
