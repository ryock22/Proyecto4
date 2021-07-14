# Proyecto4

Estudiante: Ricardo Yock

Profesor: Fabian Abarca
Curso: Modelos Probabilísticos de Señales y Sistemas
Periodo: I Semestre 2021

Asignaciones:

4.1 Modulación 16-QAM

Para la parte de modulación se separaron las portadoras por componentes seno y coseno. Y se programó 
la lógica necesaria para que el for que recorre los bits diera saltos de 4 espacios y que entre cada
salto, se evaluaran los valores de los bits i+0, i+1, i+2, i+3.
Y así sucesivamente con el resto del código. Sin embargo los resultados no fueron los esperados ya
que la imagen recuperada al final tenía muchos errores.

4.2 Estacionaridad y ergodicidad

Se sabe que una señal aleatoria es estacionaria si su comportamiento estadístico no cambia a lo largo 
del tiempo y como se puede observa en el la grafíca de "Realizaciones del proceso aleatorio" se observa
precisamente este comportamiento.

Por otro lado, se puede decir que la señal también presenta ergodicidad, porque al realizar multiples
simulaciones el comportamiento estadístico mantiene un mismo comportamiento, es decir los promedios y 
su comporamiento se mantienen temporalmente.

4.3 Densidad espectral de potencia

La densidad espectral de potencia se puede observar en la última figura la uál nos muestra la
distribución de la potencia de dicha señal sobre las distintas frecuencias que la conforman, esto nos 
permite visualizar de una forma sencilla la como los moduladores juegan con la potencia y la frecuencia
para trasmitir las señales sin que se "mezclen" o que al menos lo hagan lo menos posible.
