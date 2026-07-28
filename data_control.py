'''
data_control.py
Funciones para la lectura de la señal de un canal de ADC

V1.0 03/2026

                Axel Britos
Instituto Nacional de Tecnologia Industrial

~~~Notas de version~~~
V1.0 03/2026. Version inicial

~~~Pendiente~~~
-Agregar soporte para archivos binarios. 

'''

import numpy as np
import csv
import os
from datetime import datetime
'''
    Lee las muestras de una señal desde un archivo CSV y la devuelve como un arreglo de NumPy.
    Solamente tiene en cuenta la primer columna del archivo, ignora todo el resto.
    Soporta simbolo decimal con punto o coma. 

    Parámetros
    ----------
    file_path : str
        Ruta del archivo CSV a leer.
    max_rows : int o None, opcional
        Cantidad máxima de filas a leer. Si es None, se lee el archivo completo.
    delimiter : str, opcional
        Separador de columnas utilizado en el archivo (por defecto espacio).

    Devuelve
    --------
    numpy.ndarray
        Arreglo con los valores de la señal leídos de la primera columna del archivo.
'''
def read_from_csv(file_path, max_rows=None, delimiter = ' '):
    stat = os.stat(file_path)

    size_bytes = stat.st_size
    if size_bytes < 1024:
        size_str = f"{size_bytes} B"
    elif size_bytes < 1024**2:
        size_str = f"{size_bytes/1024:.2f} KB"
    else:
        size_str = f"{size_bytes/1024**2:.2f} MB"

    creation_time = datetime.fromtimestamp(stat.st_ctime)
    modified_time = datetime.fromtimestamp(stat.st_mtime)

    print(f"Peso             = {size_str}")
    print(f"Fecha creación   = {creation_time.strftime('%d/%m/%Y  %H:%M:%S')}")
    print(f"Última edición   = {modified_time.strftime('%d/%m/%Y  %H:%M:%S')}")

    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        signal = []
        for i, row in enumerate(reader):
            if max_rows is not None and i >= max_rows:
                break
            signal.append(float(row[0].replace(',', '.')))
    return np.array(signal)

    
#TODO implementar esto
#def read_from_raw():

def divide(signal, segment_size):
    '''
    Divide la señal para obtener la tension en cada resistor medido al multiplexar el canal. 

    El sistema inicia en configuracion Forward. Con Rx arriba y Rn abajo. La secuencia de medicion es:
    1- Forward: Tension en Rx (VRxF)
    2- Forward: Tension en Rn (VRnF)
    3- Reverse: Tension en Rn (VRnR)
    4- Reverse: Tension en Rx (VRxR)

    Se divide el vector de entrada signal en vectores de tamaño segment_size, almacenados en el array Signal_parts
    el orden el array es: VRxF, VRnF, VRnR, VRxR... y asi sucesivamente:
    VRxF = signal_parts[0], signal_parts[4], signal_parts[8]...
    VRnF = signal_parts[1], signal_parts[5], signal_parts[9]...
    VRnR = signal_parts[2], signal_parts[6], signal_parts[10]...
    VRxR = signal_parts[3], signal_parts[7], signal_parts[11]...
    
    '''
    signal_parts = [signal[i:i+segment_size] for i in range(0, len(signal), segment_size)]

    VrxF = signal_parts[0::4]
    VRnF = signal_parts[1::4]
    VRnR = signal_parts[2::4]
    VRxR = signal_parts[3::4]

    #El sistema invierte la tension al medir VRn para poder mantener la impedancia parasita constante.
    #Se invierten los valores para mantener los signos.
    VRnF = [ -1 * seg for seg in VRnF]
    VRnR = [ -1 * seg for seg in VRnR]

    return VrxF, VRnF, VRnR, VRxR


def clear_transient(signal, discard, length):
    """
    Elimina el transitorio inicial de cada segmento de señal.

    Parámetros
    ----------
    signal : list
        Lista de vectores que contienen los segmentos de señal obtenidos
        con la función divide.
    discard : int
        Cantidad de muestras iniciales a descartar en cada segmento.
    length : int
        Cantidad de muestras a conservar luego del descarte.

    Devuelve
    --------
    list
        Lista de vectores donde cada uno contiene únicamente las muestras
        seleccionadas luego de eliminar el transitorio.
    """

    cleaned_signal = []

    for segment in signal:
        cleaned_segment = segment[discard:discard + length]
        cleaned_signal.append(cleaned_segment)

    return cleaned_signal


