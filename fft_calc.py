'''
fft.py
V1.1 03/2026

                Axel Britos
Instituto Nacional de Tecnologia Industrial

~~~Notas de version~~~
V1.0 03/2026. Version inicial
V1.1 03/2026. Modifico calculate_fft() para que solo devuelva los resultados correspondientes a la mitad positiva de frecuenias.

'''
import pyfftw
import numpy as np

def calculate_fft(signal, fs):
    """
    Calcula la DFT de una señal utilizando pyFFTW y devuelve únicamente las frecuencias positivas.

    Parámetros
    ----------
    signal : np.array
        Señal en el dominio del tiempo.
    fs : float
        Frecuencia de muestreo en Hz.

    Devuelve
    ----------
    dict
        dft : numpy.ndarray (complex128)
            Espectro complejo de la señal (DFT) en formato de un solo lado (single-sided).
            Incluye la componente de continua (DC).

            La DFT está normalizada por la cantidad de muestras (1/N) y escalada por un factor 2
            en todas las componentes excepto DC, de modo que la magnitud
            |X[k]| representa directamente la amplitud de cada componente senoidal de la señal.

        freqs : numpy.ndarray
            Vector de frecuencias en Hz correspondiente a cada bin de la DFT.

        fundamental_freq_index : int
            indice de la componente de mayor magnitud (frecuencia fundamental estimada).
    """

    N = len(signal)
    input_array = pyfftw.empty_aligned(N, dtype='complex128')
    output_array = pyfftw.empty_aligned(N, dtype='complex128')
    
    input_array[:] = signal + 0j
    fft_object = pyfftw.FFTW(input_array, output_array, flags=['FFTW_ESTIMATE'])
    fft_object()
    
    #Solo me quedo con los valores de frecuencias positivas, incluyendo DC.
    half = N//2 + 1
    X = output_array[:half]

    X=X/N #Normalizacion de la FFT. 
    X[1:-1] *= 2 #Multiplico por 2 para conservar la energia total "perdida" al tomar solo la parte positiva. excepto para DC

    freqs = np.fft.fftfreq(N, d=1/fs)[:half]
    
    fundamental_freq_index = np.argmax(np.abs(X)) #La fundamental es la frecuencia con mayor magnitud.

    return {    
    "dft": X,
    "freqs": freqs,
    "fundamental_freq_index": fundamental_freq_index
    }