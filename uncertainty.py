'''
uncertainty.py
V1.1 03/2026

                Axel Britos
Instituto Nacional de Tecnologia Industrial

~~~Notas de version~~~
V1.0 03/2026. Version inicial
V1.1 03/2026. 
    *Corrijo calculo de propagacion en FFT en dft_uncertainty. 
    *agrego funcion para propagar la incertidumbre solo en frecuencia fundamental. 
V1.2 04/2026. 
    *agrego covarianza entre parte real e imaginaria en la propagacion de incertidumbre en la DFT.
    *agrego calibracion en linealidad para rango 1Vpk
V1.3 07/2026. 
    *corrijo la covarianza.
'''

import numpy as np
import sympy as sp


import matplotlib.pyplot as plt

# ---------------------------------------------
# Parámetros del ADC. ai de PXI4461
# ---------------------------------------------

n_bits = 24
V_range = 1 #Rango de entrada ADC. De -V_range a +V_range.

#---Ruido--- εn
U_noise = 2.8e-6 #Dato tomado de datasheet de PXI4461, deberia calibrarse. 
#sobre multiplicar por 6.6: Eso convertiria el valor RMS del ruido en un valor pico, me parece que no es correcto. 
#En el libro Signals and Systems de oppenheim se define que el valor rms de una señal con media seria es el desvio (sigma)

#---Cuantizacion--- εq
Q = (2 * V_range) / (2**n_bits)
U_Q = Q / np.sqrt(12) #Se divide por raiz de 12 porque se considera que la distribucion es uniforme. Para esta distribucion, el desvio es (b-a)/sqrt(12), siendo a y b los limites de la distribucion. En este caso, a=-Q/2 y b=Q/2, entonces b-a=Q.

#---Linealidad--- εlin
lin = 170e-6 #Valor obtenido de calibracion de linealidad.
U_lin = lin / np.sqrt(3) #Lo voy a sobreestimar haciendolo constante para todo el rango de medicion con distribucion rectangular..

#---Offset--- O
O = 0.0007 #Valor obtenido desde el datasheet de PXI4461. 
U_Offset = O / np.sqrt(3) #Considero que el offset es constante para todo el rango de medicion, y que tiene una distribucion uniforme.

#---Ganancia--- G
G = 0.999763605 #Maximo corrimiento al valor nominal determinado desde calibracion 
U_G = (1-G) / np.sqrt(3) #Considero que el offset es constante para todo el rango de medicion, y que tiene una distribucion uniforme.

#---Ruido 1/f--- ε1f
U_1f = 3.5e-6 #Valor obtenido del analisis mediante desvio de Allan


def pxi4461_ai_uncertainty(signal,  u_O = U_Offset, u_N = U_noise):
    '''
    Calcula la incertidumbre total de cada muestra adquirida con el módulo PXI-4461,
    considerando las contribuciones de cuantización, ruido, linealidad, offset y
    error de ganancia.

    MODELO DEL ADC - PXI4461
    x[n] = G·Vin[n] + O + εlin + εq + εn

    n = muestra
    G = ganancia [V/V]
    O = offset [V]
    εq = error de cuantificación (uniforme en [-LSB/2, LSB/2]) 
    εn = ruido [Vrms] modelado como gaussiano
    εlin = error de linealidad 

    Parámetros
    ----------
    signal : numpy.ndarray
        Vector que contiene las muestras de la señal adquirida.
    u_O : float, opcional
        Incertidumbre asociada al offset. Si el offset se calibra, corresponde a la
        incertidumbre tipo A obtenida a partir de las muestras utilizadas en la
        calibración.
    u_N : float, opcional
        Incertidumbre asociada al ruido. Si el ruido se calibra, corresponde a su
        valor RMS.

    Devuelve
    --------
    contributions : dict
        Diccionario con la contribución relativa de cada fuente de incertidumbre
        (cuantización, ruido, linealidad, offset y ganancia) respecto a la varianza total.
    signal_uncertainty : numpy.ndarray
        Vector con la incertidumbre total asociada a cada muestra de la señal.

    Notas
    -----
    Se asume independencia entre las fuentes de incertidumbre.
    '''

    #Si se envia 0 como argumento para la incertidumbre del offset o del ruido se toman los valores por defecto.  
    if u_O == 0:
        u_O = U_Offset
    if u_N == 0:
        u_N = U_noise

    #Incertidumbre total por muestra.
    signal_uncertainty = np.sqrt(
    U_Q**2 + #Cuantizacion
    u_N**2 + #Ruido
    U_lin**2 + #Linealidad
    u_O**2 + #Offset
    U_1f**2 + #Ruido 1_f
    (np.abs(signal) * U_G)**2) #Ganancia
    #TODO: Asumo que hay independencia entre las muestras. Deberia ver la correlacion para un analisis mas completo.

    #Contribucion porcentual de incertidumbre de cada termino.
    sigma_terms = {
        "Cuantizacion": U_Q,
        "Ruido": u_N,
        "Linealidad": U_lin,
        "Offset": u_O,
        "Ganancia": np.mean(np.abs(signal)) * U_G,
    }

    # ---- contribuciones porcentuales ----
    #NO MIRAR
    sigma_total_sq = sum(v**2 for v in sigma_terms.values())    # Asumo que hay independencia entre los terminos, entonces sumo las varianzas (sigma^2) para obtener la varianza total, y luego saco raiz para obtener el sigma total.
    contributions = {k: (v**2)/sigma_total_sq for k,v in sigma_terms.items()}
    
    '''
    El ε1f[n] está altamente correlacionado entre muestras consecutivas. A fs = 204.8 kS/s, dos muestras separadas por ~5 µs tienen prácticamente el mismo valor de ε1f. No son realizaciones independientes.
    Propagar ε1f directamente en el dominio de la frecuencia, donde su efecto es separable por bin. Ahí U_1f deja de ser un término por muestra y pasa a ser un término por bin de la DFT.
    '''
    return contributions, signal_uncertainty

def dft_uncertainty(u_signal, X):
    '''
    Propaga la incertidumbre de las muestras hacia el dominio frecuencial
    y calcula incertidumbre en magnitud y fase, incluyendo la covarianza
    entre parte real e imaginaria de cada bin.

    Parámetros
    ----------
    u_signal : np.array
        Incertidumbre estándar de cada muestra de la señal.
    X : np.array
        FFT de la señal (solo bins positivos).

    Devuelve
    ----------
    dict
        u_Xre  : incertidumbre parte real
        u_Xim  : incertidumbre parte imaginaria
        cov_RI : covarianza entre parte real e imaginaria
        u_A    : incertidumbre magnitud
        u_phi  : incertidumbre fase
    '''

    N = len(u_signal)
    freq_index = np.arange(len(X))

    u_Xre  = np.zeros(len(X))
    u_Xim  = np.zeros(len(X))
    cov_RI = np.zeros(len(X))
    u_A    = np.zeros(len(X))
    u_phi  = np.zeros(len(X))

    n   = np.arange(N)

    for i, k in enumerate(freq_index):

        omega_k = 2 * np.pi * k / N

        cos_vec = np.cos(omega_k * n)
        sin_vec = np.sin(omega_k * n)

        #coeficientes de sensiblidad. sale de derivar la formula de DFT respecto a cada muestra. 
        if k == 0 or (N % 2 == 0 and k == N//2):
            scale = 1 / N**2
        else:
            scale = 4 / N**2

        # --------------------------------------------------
        # Incertidumbre en Re e Im
        # --------------------------------------------------
        u_Xre[i] = np.sqrt(scale * np.sum((cos_vec**2) * u_signal**2))
        u_Xim[i] = np.sqrt(scale * np.sum((sin_vec**2) * u_signal**2))

        # --------------------------------------------------
        # Covarianza entre Re e Im
        # cada termino es dR/dx[n] * dI/dx[n] * u²[n]
        # dR/dx[n] =  (2/N)*cos(...)
        # dI/dx[n] = -(2/N)*sin(...)
        # el producto de los dos factores (2/N) es scale sin el cuadrado.
        # el signo negativo viene de dI/dx[n]
        # --------------------------------------------------
        cov_RI[i] = -scale * np.sum(cos_vec * sin_vec * u_signal**2)

        # --------------------------------------------------
        # Magnitud y fase
        # --------------------------------------------------
        Re = np.real(X[i])
        Im = np.imag(X[i])
        A  = np.sqrt(Re**2 + Im**2)

        if A < 1e-12:
            u_A[i]   = np.nan
            u_phi[i] = np.nan
            continue

        # Para propagar la incertidumbre en amplitud uso las derivadas parciales de la expresion de amplitud respecto a Re e Im.
        # dA/dRe = Re / raiz(Re^2 + Im^2) y dA/dIm = Im / raiz(Re^2 + Im^2) . Los denominadores son iguales a A.
        dA_dRe   =  Re / A
        dA_dIm   =  Im / A

        u_A[i] = np.sqrt(
            (dA_dRe**2)   * u_Xre[i]**2 +
            (dA_dIm**2)   * u_Xim[i]**2 +
            2 * dA_dRe * dA_dIm * cov_RI[i]
        )

        #La fase es atan(Im/Re). Haciendo las derivadas parciales de esa ecuacion respecto a Im Y Re se obtiene:
        dphi_dRe = -Im / (Re**2 + Im**2)
        dphi_dIm =  Re / (Re**2 + Im**2)

        u_phi[i] = np.sqrt(
            (dphi_dRe**2) * u_Xre[i]**2 +
            (dphi_dIm**2) * u_Xim[i]**2 +
            2 * dphi_dRe * dphi_dIm * cov_RI[i]
        )

    return {
        "u_Xre":  u_Xre,
        "u_Xim":  u_Xim,
        "cov_RI": cov_RI,
        "u_A":    u_A,
        "u_phi":  u_phi,
    }

def dft_uncertainty_fundamental(u_signal, X, k_fund):
    '''
    Propaga la incertidumbre de las muestras hacia la frecuencia fundamental,
    incluyendo covarianza entre parte real e imaginaria.

    Parámetros
    ----------
    u_signal : np.array
        Incertidumbre estándar de cada muestra de la señal.
    X : np.array
        FFT de la señal (solo bins positivos).
    k_fund : int
        Índice del bin de la frecuencia fundamental.

    Devuelve
    ----------
    dict
        u_Xre  : incertidumbre parte real
        u_Xim  : incertidumbre parte imaginaria
        cov_RI : covarianza entre parte real e imaginaria
        u_A    : incertidumbre magnitud
        u_phi  : incertidumbre fase
    '''

    N   = len(u_signal)
    n   = np.arange(N)
    u2  = u_signal**2

    omega_k = 2 * np.pi * k_fund / N
    cos_vec = np.cos(omega_k * n)
    sin_vec = np.sin(omega_k * n)

    #coeficientes de sensiblidad. sale de derivar la formula de DFT respecto a las muestras.
    #En realidad es 1/N y 2/N pero ya los dejo elevado al cuadrado porque se usan asi para el calculo de la incertidumbre.
    if k_fund == 0 or (N % 2 == 0 and k_fund == N//2):
        scale = 1 / N**2
    else:
        scale = 4 / N**2

    # --------------------------------------------------
    # Incertidumbre en Re e Im
    # --------------------------------------------------
    u_Xre  = np.sqrt(scale * np.sum((cos_vec**2) * u2))
    u_Xim  = np.sqrt(scale * np.sum((sin_vec**2) * u2))

    # --------------------------------------------------
    # Covarianza entre Re e Im
    # --------------------------------------------------
    cov_RI = -np.sum(scale *cos_vec *sin_vec * u2)

    # --------------------------------------------------
    # Magnitud y fase
    # --------------------------------------------------
    Re = np.real(X[k_fund])
    Im = np.imag(X[k_fund])
    A  = np.sqrt(Re**2 + Im**2)

    if A < 1e-12:
        return {
            "u_Xre":  u_Xre,
            "u_Xim":  u_Xim,
            "cov_RI": cov_RI,
            "u_A":    np.nan,
            "u_phi":  np.nan,
        }

    # Para propagar la incertidumbre en amplitud uso las derivadas parciales de la expresion de amplitud respecto a Re e Im.
    # dA/dRe = Re / raiz(Re^2 + Im^2) y dA/dIm = Im / raiz(Re^2 + Im^2) . Los denominadores son iguales a A.
    dA_dRe   =  Re / A
    dA_dIm   =  Im / A

    u_A = np.sqrt(
        (dA_dRe**2)   * u_Xre**2 +
        (dA_dIm**2)   * u_Xim**2 +
        2 * dA_dRe * dA_dIm * cov_RI
    )
    #La fase es atan(Im/Re). Haciendo las derivadas parciales de esa ecuacion respecto a Im Y Re se obtiene:
    dphi_dRe = -Im / (Re**2 + Im**2) #derivada de la fase respecto a Re
    dphi_dIm =  Re / (Re**2 + Im**2) #derivada de la fase respecto a Im


    u_phi = np.sqrt(
        (dphi_dRe**2) * u_Xre**2 +
        (dphi_dIm**2) * u_Xim**2 +
        2 * dphi_dRe * dphi_dIm * cov_RI
    )

    return {
        "u_Xre":  u_Xre,
        "u_Xim":  u_Xim,
        "cov_RI": cov_RI,
        "u_A":    u_A,
        "u_phi":  u_phi,
    }
def correction_uncertainty(voltages, uncertainties):
    """
    Calcula la propagación de incertidumbre del factor de corrección complejo k,
    a partir de cuatro tensiones complejas medidas (VRxF, VRnF, VRnR, VRxR) y
    sus incertidumbres asociadas en parte real e imaginaria.

    Parameters
    ----------
    voltages : iterable
        Iterable de longitud N, donde cada elemento contiene 4 tensiones complejas:
            (VRxF, VRnF, VRnR, VRxR)

    uncertainties : iterable
        Iterable de longitud N, donde cada elemento contiene las incertidumbres
        asociadas a cada tensión compleja, expresadas como:
            (
                (u_x1, u_y1, cov1),
                (u_x2, u_y2, cov2),
                (u_x3, u_y3, cov3),
                (u_x4, u_y4, cov4)
            )

    Returns
    -------
    dict
        Diccionario con las siguientes claves:
            - "u_k_mag"
            - "u_k_phase"
    """

    #--------------------------------------------------
    # 1. Definicion de variables simbolicas
    #--------------------------------------------------
    x1, y1 = sp.symbols('x1 y1', real=True)
    x2, y2 = sp.symbols('x2 y2', real=True)
    x3, y3 = sp.symbols('x3 y3', real=True)
    x4, y4 = sp.symbols('x4 y4', real=True)

    VRxF = x1 + sp.I*y1
    VRnF = x2 + sp.I*y2
    VRnR = x3 + sp.I*y3
    VRxR = x4 + sp.I*y4

    #--------------------------------------------------
    # 2. Fórmula de corrección
    #--------------------------------------------------
    A = VRxF / VRnF
    B = VRnR / VRxR
    k = (A + 1) / (B + 1)

    k_re = sp.re(k)
    k_im = sp.im(k)

    vars = (x1, y1, x2, y2, x3, y3, x4, y4)

    dk_re = [sp.diff(k_re, v) for v in vars]
    dk_im = [sp.diff(k_im, v) for v in vars]

    f_re = [sp.lambdify(vars, expr, 'numpy') for expr in dk_re]
    f_im = [sp.lambdify(vars, expr, 'numpy') for expr in dk_im]

    #--------------------------------------------------
    # 3. Inicialización
    #--------------------------------------------------
    u_k_re = []
    u_k_im = []
    cov_k_re_im = []

    u_k_mag = []
    u_k_phase = []

    #--------------------------------------------------
    # 4. Loop para segmento
    #--------------------------------------------------
    for v, u in zip(voltages, uncertainties):

        VRxF_val, VRnF_val, VRnR_val, VRxR_val = v

        (u_x1, u_y1, cov1), \
        (u_x2, u_y2, cov2), \
        (u_x3, u_y3, cov3), \
        (u_x4, u_y4, cov4) = u

        values = (
            np.real(VRxF_val), np.imag(VRxF_val),
            np.real(VRnF_val), np.imag(VRnF_val),
            np.real(VRnR_val), np.imag(VRnR_val),
            np.real(VRxR_val), np.imag(VRxR_val)
        )

        #--------------------------------------------------
        # Covarianza de las variables de entrada
        #--------------------------------------------------
        Ux = np.array([

            [u_x1**2, cov1,      0,       0,       0,       0,       0,       0],
            [cov1,    u_y1**2,   0,       0,       0,       0,       0,       0],

            [0,       0,         u_x2**2, cov2,    0,       0,       0,       0],
            [0,       0,         cov2,    u_y2**2, 0,       0,       0,       0],

            [0,       0,         0,       0,       u_x3**2, cov3,    0,       0],
            [0,       0,         0,       0,       cov3,    u_y3**2, 0,       0],

            [0,       0,         0,       0,       0,       0,       u_x4**2, cov4],
            [0,       0,         0,       0,       0,       0,       cov4,    u_y4**2]

        ], dtype=float)

        #--------------------------------------------------
        # Jacobiano
        #--------------------------------------------------
        J = np.array([
            [f(*values) for f in f_re],
            [f(*values) for f in f_im]
        ], dtype=float)

        #--------------------------------------------------
        # Propagación según GUM
        #--------------------------------------------------
        Uk = J @ Ux @ J.T

        var_re = Uk[0, 0]
        var_im = Uk[1, 1]
        cov_re_im = Uk[0, 1]

        uk_re = np.sqrt(var_re)
        uk_im = np.sqrt(var_im)

        u_k_re.append(uk_re)
        u_k_im.append(uk_im)
        cov_k_re_im.append(cov_re_im)

        #--------------------------------------------------
        # Cálculo de k
        #--------------------------------------------------
        k_val = (VRxF_val / VRnF_val + 1) / (VRnR_val / VRxR_val + 1)

        Re = np.real(k_val)
        Im = np.imag(k_val)

        A = np.sqrt(Re**2 + Im**2)

        if A < 1e-12:
            u_k_mag.append(np.nan)
            u_k_phase.append(np.nan)
            continue

        #--------------------------------------------------
        # Propagación de incertidumbre a coordenadas polares
        #--------------------------------------------------
        dA_dRe = Re / A
        dA_dIm = Im / A

        dphi_dRe = -Im / (Re**2 + Im**2)
        dphi_dIm = Re / (Re**2 + Im**2)

        uA = np.sqrt(
            (dA_dRe**2) * var_re +
            (dA_dIm**2) * var_im +
            2 * dA_dRe * dA_dIm * cov_re_im
        )

        uPhi = np.sqrt(
            (dphi_dRe**2) * var_re +
            (dphi_dIm**2) * var_im +
            2 * dphi_dRe * dphi_dIm * cov_re_im
        )

        u_k_mag.append(uA)
        u_k_phase.append(uPhi)

    #--------------------------------------------------
    # 5. Salida
    #--------------------------------------------------
    return {
        #"u_k_re": np.array(u_k_re),
        #"u_k_im": np.array(u_k_im),
        #"cov_k_re_im": np.array(cov_k_re_im),
        "u_k_mag": np.array(u_k_mag),
        "u_k_phase": np.array(u_k_phase)
    }

def calc_Rx(k_values, k_unc, Rn_mag, Rn_phase, u_Rn_mag, u_Rn_phase):
    """
    Calcula el valor complejo de un resistor desconocido Rx a partir de un resistor
    de referencia Rn y un factor complejo k, propagando la incertidumbre en magnitud
    y fase.

    Modelo utilizado:
        |Rx| = |Rn| * |k|
        φ_Rx = φ_Rn + φ_k

    Parámetros
    ----------
    k_values : array-like of complex
        Valores complejos del factor k para cada medición.

    k_unc : dict
        Diccionario con las incertidumbres de k:
            - "u_k_mag": array-like, incertidumbre de la magnitud de k
            - "u_k_phase": array-like, incertidumbre de la fase de k [rad]

    Rn_mag : float
        Magnitud del resistor de referencia Rn.

    Rn_phase : float
        Fase del resistor de referencia Rn [rad].

    u_Rn_mag : float
        Incertidumbre de la magnitud de Rn.

    u_Rn_phase : float
        Incertidumbre de la fase de Rn [rad].

    Retorna
    -------
    dict
        Diccionario con los resultados:
            - "Rx_mag": ndarray, magnitud de Rx
            - "Rx_phase": ndarray, fase de Rx [rad]
            - "u_Rx_mag": ndarray, incertidumbre de la magnitud
            - "u_Rx_phase": ndarray, incertidumbre de la fase
    """
    Rx_mag = []
    Rx_phase = []

    u_Rx_mag = []
    u_Rx_phase = []

    for k, u_k_mag, u_k_phase in zip(
        k_values,
        k_unc["u_k_mag"],
        k_unc["u_k_phase"]
    ):

        #Expresa factor k en coordenadas polares.
        k_mag = np.abs(k)
        k_phase = np.angle(k)

        #-------------------------
        # Magnitud
        #-------------------------
        Rx_m = Rn_mag * k_mag #Obtengo la magnitud de Rx

        #Hago las derivadas parciales de la expresion de Rx_m
        dRx_dRn = k_mag #Derivada de Rx_m respecto a Rn_mag 
        dRx_dk = Rn_mag#Derivada de Rx_m respecto a k_mag
        
        #Propago con la incertidumbre de Rn y de k
        u_m=np.sqrt((dRx_dRn * u_Rn_mag)**2 + (dRx_dk * u_k_mag)**2)
   
        #-------------------------
        # Fase
        #-------------------------
        Rx_p = Rn_phase + k_phase #Obtengo la fase de Rx
        #Hago las derivadas parciales de la expresion de Rx_p
        #dRx_p/dRn_phase = 1
        #dRx_p/dk_phase = 1

        #propago
        u_p = np.sqrt(
            u_Rn_phase**2 +
            u_k_phase**2
        )

        #------------------------
        # genero los arrays
        #------------------------
        Rx_mag.append(Rx_m)
        Rx_phase.append(Rx_p)
        u_Rx_mag.append(u_m)
        u_Rx_phase.append(u_p)

    return {
        "Rx_mag": np.array(Rx_mag),
        "Rx_phase": np.array(Rx_phase),
        "u_Rx_mag": np.array(u_Rx_mag),
        "u_Rx_phase": np.array(u_Rx_phase)
    }