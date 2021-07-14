# -*- coding: utf-8 -*-
"""
Julio 2021 - I Semestre

@author: Ricardo Yock
"""

# Proyecto 4 Modelos probabilísticos de señales y sistemas (IE-0405)

# Librerías utilizadas

from PIL import Image 
import numpy as np 
import matplotlib.pyplot as plt
import time 
from scipy import fft 



#**************************************************************************
# 4.1 Modulación 16-QAM
#**************************************************************************

def fuente_info(imagen):
    # Convertimos la imagen en un arreglo de pixeles utilizando la librería de nunmpy
    
    img = Image.open(imagen)
    
    return np.array(img)


def rgb_a_bit(imagen):
    '''Convertimos los pixeles de base decimal a binaria. '''
       
    x, y, z = imagen.shape                  # Dimensiones de la imagen
    n_pixeles = x * y * z                   # Número total de pixeles
    pixeles = np.reshape(imagen, n_pixeles) # Vector unidimensional de n_pixeles

    # Convertir los canales a base 2
    bits = [format(pixel,'08b') for pixel in pixeles]
    bits_Rx = np.array(list(''.join(bits)))
    
    return bits_Rx.astype(int)


''' Se separará el modulador seno y coseno por simplicidad.''' 

# Moduladora coseno

def modcoseno(bits, fc, mpp):
   
    N = len(bits) # Cantidad de bits
    Tc = 1 / fc  
    t_periodo = np.linspace(0, Tc, mpp)
    portcoseno = np.cos(2*np.pi*fc*t_periodo)

    # Inicializar la señal modulada s(t)
    t_simulacion     = np.linspace(0, N*Tc, N*mpp) 
    senal_Txcoseno   = np.zeros(t_simulacion.shape)
    moduladoracoseno = np.zeros(t_simulacion.shape)  # Señal de información
 
    # Formas de onda según los bits
    for i in range(0, len(bits), 4):
        if bits[i] == 0 and bits[i+1] == 0:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * -3
            moduladoracoseno[i*mpp : (i+1)*mpp] = 0

        elif bits[i] == 0 and bits[i+1] == 1:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * -1
            moduladoracoseno[i*mpp : (i+1)*mpp] = 0

        elif bits[i] == 1 and bits[i+1] == 1:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * 1
            moduladoracoseno[i*mpp : (i+1)*mpp] = 1

        elif bits[i] == 1 and bits[i+1] == 0:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * 3
            moduladoracoseno[i*mpp : (i+1)*mpp] = 1
        else:
            pass
            
    # Potencia promedio de la señal modulada
    Pmcoseno = (1 / (N*Tc)) * np.trapz(pow(senal_Txcoseno, 2), t_simulacion)

    return senal_Txcoseno, Pmcoseno, portcoseno, moduladoracoseno       
            


# Moduladora seno

def modseno(bits, fc, mpp):
   
    N = len(bits) # Cantidad de bits
    Tc = 1 / fc  
    t_periodo = np.linspace(0, Tc, mpp)
    portseno = np.cos(2*np.pi*fc*t_periodo)

    # Inicializar la señal modulada s(t)
    t_simulacion   = np.linspace(0, N*Tc, N*mpp) 
    senal_Txseno   = np.zeros(t_simulacion.shape)
    moduladoraseno = np.zeros(t_simulacion.shape)  # señal de información

    # Formas de onda según los bits
    for i in range(0, len(bits), 4):
        if bits[i+2] == 0 and bits[i+3] == 0:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * 3
            moduladoracoseno[i*mpp : (i+1)*mpp] = 1

        elif bits[i+2] == 0 and bits[i+3] == 1:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * 1
            moduladoracoseno[i*mpp : (i+1)*mpp] = 1

        elif bits[i+2] == 1 and bits[i+3] == 1:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * -1
            moduladoracoseno[i*mpp : (i+1)*mpp] = 0

        elif bits[i+2] == 1 and bits[i+3] == 0:
            senal_Txcoseno[i*mpp : (i+1)*mpp] = portcoseno * -3
            moduladoracoseno[i*mpp : (i+1)*mpp] = 0
        else:
            pass

    # Potencia promedio de la señal modulada
    Pmseno = (1 / (N*Tc)) * np.trapz(pow(senal_Txseno, 2), t_simulacion)
    
    return senal_Txseno, Pmseno, portseno, moduladoraseno 



def canal_ruidoso(signal_Txx, Pm, SNR):
    
    # Un bloque que simula ruido AWGN.  
    
    Pn = Pm / pow(10, SNR/10)                                  # Potencia del ruido generado por el canal
    ruido = np.random.normal(0, np.sqrt(Pn), signal_Txx.shape) # Ruido AWGN
    senal_Rx = signal_Txx + ruido                              # Señal distorsionada por el canal ruidoso
   
    return senal_Rx


def demodulador(senal_Rx, portcoseno, portseno, mpp):
    
    # El demodulador decodifica la señal mediante la detección de energía.
     
    portadora = portcoseno + portseno # Suma de las dos señales portadoras ortogonales.
    M = len(senal_Rx) # Cantidad de muestras en senal_Rx
    N = int(M / mpp) # Cantidad de bits en transmisión
    bits_Rx = np.zeros(N) # Vector para bits obtenidos por la demodulación
    senal_demodulada = np.zeros(M) # Vector para la señal demodulada
    Es = np.sum(portadora**2) # Energía de un período de la portadora

    # Demodulación
    for i in range(N):
        # Producto interno de dos funciones
        producto = senal_Rx[i*mpp : (i+1)*mpp] * portcoseno + senal_Rx[i*mpp : (i+1)*mpp] * portseno
        senal_demodulada[i*mpp : (i+1)*mpp] = producto
        Ep = np.sum(producto) 

        # Criterio de decisión por detección de energía
        if Ep > Es*0.8:
            bits_Rx[i] = 1
        else:
            bits_Rx[i] = 0

    return bits_Rx.astype(int), senal_demodulada




def bits_a_rgb(bits_Rx, dimensiones):
    # Función que decodifica los bits recuperados en el proceso de demodulación

    '''
    :param: Un vector de bits 1 x k 
    :param dimensiones: Tupla con dimensiones de la img.
    :return: Un array con los pixeles reconstruidos
    '''
    # Cantidad de bits
    N = len(bits_Rx)

    # Se reconstruyen los canales RGB
    bits = np.split(bits_Rx, N / 8)

    # Se decofican los canales:
    canales = [int(''.join(map(str, canal)), 2) for canal in bits]
    pixeles = np.reshape(canales, dimensiones)

    return pixeles.astype(np.uint8)


''' Entrega de información '''


# Parámetros
fc = 5000  # frecuencia de la portadora
mpp = 40   # muestras por periodo de la portadora
SNR = 0    # relación señal-a-ruido del canal

# Iniciar medición del tiempo de simulación
inicio = time.time()

# 1. Importar y convertir la imagen a trasmitir
imagen_Tx = fuente_info('arenal.jpg')
dimensiones = imagen_Tx.shape

# 2. Codificar los pixeles de la imagen
bits_Tx = rgb_a_bit(imagen_Tx)

# 3.1 Moduladora coseno 
senal_Txcoseno, Pmcoseno, portcoseno, moduladoracoseno= modcoseno(bits_Tx, fc, mpp)

# 3.2 Moduladora seno
senal_Txseno, Pmseno, portseno, moduladoraseno= modseno(bits_Tx, fc, mpp)

#definiciones 

portadora  = portcoseno + portseno
Pm         = Pmcoseno + Pmseno
signal_Txx = senal_Txcoseno + senal_Txseno
moduladora = moduladoracoseno + moduladoraseno

# 4. Se transmite la señal modulada, por un canal ruidoso
senal_Rx = canal_ruidoso(signal_Txx, Pm, SNR)

# 5. Se desmodula la señal recibida del canal
bits_Rx, senal_demodulada = demodulador(senal_Rx, portcoseno, portseno, mpp)

# 6. Se visualiza la imagen recibida 
imagen_Rx = bits_a_rgb(bits_Rx, dimensiones)
Fig = plt.figure(figsize=(10,6))

# Cálculo del tiempo de simulación
print('Duración de la simulación: ', time.time() - inicio)

# 7. Calcular número de errores
errores = sum(abs(bits_Tx - bits_Rx))
BER = errores/len(bits_Tx)
print('{} errores, para un BER de {:0.4f}.'.format(errores, BER))

# Mostrar imagen transmitida
ax = Fig.add_subplot(1, 2, 1)
imgplot = plt.imshow(imagen_Tx)
ax.set_title('Transmitido')

# Mostrar imagen recuperada
ax = Fig.add_subplot(1, 2, 2)
imgplot = plt.imshow(imagen_Rx)
ax.set_title('Recuperado')
Fig.tight_layout()

plt.imshow(imagen_Rx)

# Visualizar el cambio entre las señales
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, sharex=True, figsize=(14, 7))

# La onda cuadrada moduladora (bits de entrada)
ax1.plot(moduladora[0:600], color='r', lw=2) 
ax1.set_ylabel('$b(t)$')

# La señal modulada con 16-QAM
ax2.plot(signal_Txx[0:600], color='g', lw=2) 
ax2.set_ylabel('$s(t)$')

# La señal modulada al dejar el canal
ax3.plot(senal_Rx[0:600], color='b', lw=2) 
ax3.set_ylabel('$s(t) + n(t)$')

# La señal demodulada
ax4.plot(senal_demodulada[0:600], color='m', lw=2) 
ax4.set_ylabel('$b^{\prime}(t)$')
ax4.set_xlabel('$t$ / milisegundos')
fig.tight_layout()
plt.show()






#**************************************************************************
# 4.2. Estacionariedad y ergodicidad 
#**************************************************************************



t_final   = 0.1
T_per     = 100	# tiempo en segundos
tmuestreo = np.linspace(0, t_final, T_per)

# Inicialización del proceso aleatorio X(t) con Nr = 4 realizaciones
realizaciones = 4
x_t = np.empty((realizaciones, len(tmuestreo)))	# Nr = 4 funciones del tiempo x(t) con T puntos

# Creación de las muestras del proceso x(t) (A y Z independientes)
A_j = [-1,1]
vacos =  np.cos(2*(np.pi)*fc*tmuestreo)
vasen =  np.sin(2*(np.pi)*fc*tmuestreo)

for i in A_j:
    X = i*vacos + i*vasen 
    Y = -i*vacos + i*vasen 
    x_t[i,:] = X
    x_t[i+1,:] = Y
    plt.plot(tmuestreo, X)
    plt.plot(tmuestreo, Y)


# Promedio de las N realizaciones en cada instante (cada punto en t)
P = [np.mean(x_t[:,i]) for i in range(len(tmuestreo))]
plt.plot(tmuestreo, P, lw=6, color = '#0B5394', label = 'Promedio')

# Graficar el resultado teórico del valor esperado
E = np.mean(signal_Txx)*tmuestreo
plt.plot(tmuestreo, E, '-.', lw=4, color='#EA9999', label = 'Teórico')

# Mostrar las realizaciones, y su promedio calculado y teórico
plt.title('Realizaciones del proceso aleatorio $X(t)$')
plt.xlabel('$t$')
plt.ylabel('$x_i(t)$')
plt.show()






#**************************************************************************
# 4.3 Densidad espectral de potencia 
#**************************************************************************

fourier = fft(signal_Txx)

Nm = len(signal_Txx) # Muestras de la señal
Ns = Nm // mpp       # Número de símbolos 
Tc = 1/ fc           # Tiempo del símbolo = período de la onda portadora
Tm = Tc / mpp        # Período de muestreo
T = Ns * Tc          # Tiempo de la simulación
f = np.linspace(0.0, 1.0/(2.0*Tm), Nm//2) # Espacio de frecuencias

# Gráfica
plt.plot(f,2.0/Nm*np.power(np.abs(fourier[0:Nm//2]),2), color = '#EA9999')
plt.xlim(0,30000)
plt.grid()
plt.show()


