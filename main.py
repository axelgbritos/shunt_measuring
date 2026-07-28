'''


Medicion de Fase

Procesador

V1.2 03/2026

Axel Britos
Instituto Nacional de Tecnologia Industrial

Notas de version:
V1.1- 
    *elimino la opcion de activar o desactivar la impresion de datos. 
    *mensajes se imprimen en pantalla con consola rich
V1.2-
    *compatibilidad con version 1.1 de uncertainty.py. Hago la propagacion de incertidumbre a DFT solo para frecuencia fundamental. 
'''

import argparse
import sys
import time

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn, TextColumn

import data_control
import uncertainty
import fft_calc
from grapher import pageGraphs
import tables
import tools

#import linker

grapher = pageGraphs()
console = Console()

#--------------------------------------------------------------------
# 0 - Configurar parametros de programa.
#--------------------------------------------------------------------
GRAPH_INFO = False                 #Grafica los resultados parciales del proceso.
DEFAULT_MEASURE_PATH = "RAW/Rcal10-vs-1tanos/1000hz-1.csv"         #Direccion del archivo de medicion.
#MEASURE_PATH = "RAW/11-25-6-10vs2.csv"         #Direccion del archivo de medicion.

#---Resistor de Referencia 
DEFAULT_Rn_mag = 10.148
DEFAULT_U_Rn_mag = 100e-6

DEFAULT_Rn_pha = 1.12E-04 
DEFAULT_U_Rn_pha	=1.48E-06


fs = 204800                         #Frecuencia de muestreo [Muestras/segundo]
f_switching = 10                    #Frecuencia de cambio de configuracion [Hz]
transient_lenght = 2048*2           #Cantidad de muestras que se descartan, correspondiente al transitorio.
signal_lenght = 2048*6              #Cantidad de muestras limpias de la señal.
valid_segments = 5                  #Cantidad de segmentos que seran analizados. La variable representa los ultimos segmentos a analizar.

#---Calibraciones previas: 
# Offset. 
offset = 1.41645E-05    #Si se calibro el offset ingresar el valor medio de la señal de calibracion. Sino, poner 0
sigma_offset = 6.88605E-09;   #Si se calibro el offset ingresar la incertidumbre A de la señal de calibracion. (desv / raiz(n)). Sino, poner 0
# Ruido. 
sigma_noise = 3.11626E-06;    #Si se calibro el offset ingresar el valor RMS de la señal de calibracion menos el valor medio. Sino, poner 0

'''
si se abrio desde linker, reasignar valores de las variables...
'''

def main(
    Rn_mag=DEFAULT_Rn_mag,
    U_Rn_mag=DEFAULT_U_Rn_mag,
    Rn_pha=DEFAULT_Rn_pha,
    U_Rn_pha=DEFAULT_U_Rn_pha,
    MEASURE_PATH=DEFAULT_MEASURE_PATH,
):
    #--------IGNORAR. Interface grafica----------------
    with Progress(
        SpinnerColumn(style="cyan"),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Procesando...", total=None)
    #---------------------IGNORAR----------------------

    #--------------------------------------------------------------------
    # 1- Abrir archivo de medicion
    #--------------------------------------------------------------------
        console.rule("[bold green]Seleccion de las muestras")
        progress.update(task, description="[cyan]Leyendo archivo")
        
        signal = data_control.read_from_csv(MEASURE_PATH, max_rows=None)
        console.print(f"[bold]Archivo {MEASURE_PATH} leido. Muestras: {len(signal)}")
    #--------------------------------------------------------------------
    # 2- Cortar la señal y obtener la tension en cada resistor en cada configuracion.
    #--------------------------------------------------------------------
        progress.update(task, description="[cyan]Ordenando señales")
        
        segment_size = fs // f_switching    #Cantidad de muestras que tendra cada segmento obtenido.
        VRxF, VRnF, VRnR, VRxR = data_control.divide(signal,segment_size)

        console.print("Señales separadas.")

        signalarrays = {
            "VRxF": VRxF,
            "VRxR": VRxR,
            "VRnF": VRnF,
            "VRnR": VRnR
        }
        for name, data in signalarrays.items():
            n_vectors = len(data)
            n_samples = len(data[0]) if n_vectors > 0 else 0

            console.print(f"{name}: {n_vectors} vectores de {n_samples} muestras cada uno.")

        #Grafico de la señal de entrada. 
        grapher.signal_segments(signal, fs//f_switching, transient_lenght, signal_lenght, titulo=" ")
    #--------------------------------------------------------------------
    # 3- Limpiar transitorios
    #--------------------------------------------------------------------
        progress.update(task, description="[cyan]Removiendo Transitorio")

        VRxF = VRxF[-valid_segments:] 
        VRnF = VRnF[-valid_segments:]
        VRnR = VRnR[-valid_segments:]
        VRxR = VRxR[-valid_segments:]

        VRxF=data_control.clear_transient(VRxF,transient_lenght,signal_lenght)
        VRnF=data_control.clear_transient(VRnF,transient_lenght,signal_lenght)
        VRnR=data_control.clear_transient(VRnR,transient_lenght,signal_lenght)
        VRxR=data_control.clear_transient(VRxR,transient_lenght,signal_lenght)

        console.print("Transitorio removido.")

        for name, data in signalarrays.items():
            n_vectors = len(data)
            n_samples = len(data[0]) if n_vectors > 0 else 0

            console.print(f"{name}: {n_vectors} vectores de {n_samples} muestras cada uno.")

    #--------------------------------------------------------------------
    # 4- Compensar Offset (si corresponde) a cada set de datos.
    #--------------------------------------------------------------------
        if offset != 0: #Se corrige el offset si se ingreso como dato.
            progress.update(task, description="[cyan]Compensando Offset")
            VrxF = [segment - offset for segment in VRxF]
            VRnF = [segment - offset for segment in VRnF]
            VRnR = [segment - offset for segment in VRnR]
            VRxR = [segment - offset for segment in VRxR]

            console.print(f"Correccion de offset aplicada: {tools.format_with_prefix(offset,unit='V')}.")

    #--------------------------------------------------------------------
    # 5- Calcular incertidumbre a cada set de datos.
    #--------------------------------------------------------------------
        console.rule("[bold green]Incertidumbre - Dominio Temporal")
        progress.update(task, description="[cyan]Calculando Incertidumbre a las muestras")

        #Ver el modelo del ADC y sus incertidumbres en uncertainty.py
        U_VRxF = [uncertainty.pxi4461_ai_uncertainty(segment,u_O= sigma_offset,u_N=sigma_noise)[1] for segment in VRxF]
        U_VRxR = [uncertainty.pxi4461_ai_uncertainty(segment,u_O= sigma_offset,u_N=sigma_noise)[1] for segment in VRxR]
        U_VRnF = [uncertainty.pxi4461_ai_uncertainty(segment,u_O= sigma_offset,u_N=sigma_noise)[1] for segment in VRnF]
        U_VRnR = [uncertainty.pxi4461_ai_uncertainty(segment,u_O= sigma_offset,u_N=sigma_noise)[1] for segment in VRnR]

        console.print("\nFuentes de Incertidumbre: ")
        console.print(f"Ganancia: {tools.format_with_prefix(uncertainty.U_G)}")
        console.print(f"Ruido: {tools.format_with_prefix(sigma_noise)}")
        console.print(f"Offset: {tools.format_with_prefix(sigma_offset)}")
        console.print(f"Linealidad: {tools.format_with_prefix(uncertainty.U_lin)}")
        console.print(f"Cuantizacion: {tools.format_with_prefix(uncertainty.U_Q)}")
        console.print(f"Ruido 1/f: {tools.format_with_prefix(uncertainty.U_1f)}\n")

        console.print("Contribución porcentual de cada componente:")
        contributions, _ = uncertainty.pxi4461_ai_uncertainty(VRxF[0],u_O=sigma_offset,u_N=sigma_noise)
        for k, v in contributions.items():
            console.print(f"{k}: {100*v:.4f}%")

        #Grafico de histogramas de incertidumbre.
        def hist1(ax): ax.hist(U_VRxF[0], bins=50); ax.set_title("Incertidumbre en VRxF"); ax.set_xlabel("Valor"); ax.set_ylabel("Frecuencia")
        def hist2(ax): ax.hist(U_VRnF[0], bins=50); ax.set_title("Incertidumbre en VRnF"); ax.set_xlabel("Valor"); ax.set_ylabel("Frecuencia")
        grapher.add_subplots([[hist1, hist2]], titulo="")

    #--------------------------------------------------------------------
    # 6- Calcular DFT de cada set de datos. 
    #--------------------------------------------------------------------
        progress.update(task, description="[cyan]Calculando DFT.")

        VRxF_dft = [fft_calc.calculate_fft(segment, fs=fs) for segment in VRxF]
        VRxR_dft = [fft_calc.calculate_fft(segment, fs=fs) for segment in VRxR]
        VRnF_dft = [fft_calc.calculate_fft(segment, fs=fs) for segment in VRnF]
        VRnR_dft = [fft_calc.calculate_fft(segment, fs=fs) for segment in VRnR]

    #--------------------------------------------------------------------
    # 7- Propagar incertidumbre a DFT
    #--------------------------------------------------------------------   
        progress.update(task, description="[cyan]Propagando incertidumbre a la frecuencia fundamental de DFT.")

        U_VRxF_dft = [uncertainty.dft_uncertainty_fundamental(u_signal=u, X=d["dft"],k_fund=d["fundamental_freq_index"]) for u, d in zip(U_VRxF, VRxF_dft)]
        U_VRxR_dft = [uncertainty.dft_uncertainty_fundamental(u_signal=u, X=d["dft"],k_fund=d["fundamental_freq_index"]) for u, d in zip(U_VRxR, VRxR_dft)]
        U_VRnF_dft = [uncertainty.dft_uncertainty_fundamental(u_signal=u, X=d["dft"],k_fund=d["fundamental_freq_index"]) for u, d in zip(U_VRnF, VRnF_dft)]
        U_VRnR_dft = [uncertainty.dft_uncertainty_fundamental(u_signal=u, X=d["dft"],k_fund=d["fundamental_freq_index"]) for u, d in zip(U_VRnR, VRnR_dft)]

        fund_freq = VRxF_dft[0]["freqs"][VRxF_dft[0]["fundamental_freq_index"]]
        print(f"Frecuencia fundamental: {tools.format_with_prefix(fund_freq, 'Hz')}")

        console.rule("[bold green]Tensiones")
        console.print("[italic white]Las desviaciones son solamente sobre los valores mostrados en la tabla. En este paso no toma como incertidumbre.")
        tables.print_fft_table(console, VRxF_dft, U_VRxF_dft, title="Tension en Rx - Forward")
        tables.print_fft_table(console, VRxR_dft, U_VRxR_dft, title="Tension en Rx - Reverse")
        tables.print_fft_table(console, VRnF_dft, U_VRnF_dft, title="Tension en Rn - Forward")
        tables.print_fft_table(console, VRnR_dft, U_VRnR_dft, title="Tension en Rn - Reverse")

    #--------------------------------------------------------------------
    # 8- Aplicar algoritmo de correccion
    #--------------------------------------------------------------------
        console.rule("[bold green]Correccion")
        progress.update(task, description="[cyan]Aplicando Correccion")

        k_values = [] #lista para los valores de los valores con correccion aplicada.

        for vrxf, vrnf, vrnr, vrxr in zip(VRxF_dft, VRnF_dft, VRnR_dft, VRxR_dft):

            idx = vrxf["fundamental_freq_index"]  # mismo índice para todos

            VRxF_val = vrxf["dft"][idx]
            VRnF_val = vrnf["dft"][idx]
            VRnR_val = vrnr["dft"][idx]
            VRxR_val = vrxr["dft"][idx]

            A = VRxF_val / VRnF_val
            B = VRnR_val / VRxR_val

            k1 = A + 1
            k2 = B + 1

            k = k1 / k2

            k_values.append(k)
    #--------------------------------------------------------------------
    # 9- Propagar incertidumbre a formula de correccion.
    #--------------------------------------------------------------------
        progress.update(task, description="[cyan]Propagando incertidumbre a la formula de correccion")

        # Generacion de la lista de tensiones complejas.
        voltages = []
        for vrxf, vrnf, vrnr, vrxr in zip(VRxF_dft, VRnF_dft, VRnR_dft, VRxR_dft):
            idx = vrxf["fundamental_freq_index"]

            voltages.append((
                vrxf["dft"][idx],
                vrnf["dft"][idx],
                vrnr["dft"][idx],
                vrxr["dft"][idx]
            ))

        # Generacion de la lista de incertidumbres
        uncertainties = []
        for u_vrxf, u_vrnf, u_vrnr, u_vrxr in zip(
            U_VRxF_dft,
            U_VRnF_dft,
            U_VRnR_dft,
            U_VRxR_dft
        ):

            uncertainties.append((
                (u_vrxf["u_Xre"], u_vrxf["u_Xim"], u_vrxf["cov_RI"]),
                (u_vrnf["u_Xre"], u_vrnf["u_Xim"], u_vrnf["cov_RI"]),
                (u_vrnr["u_Xre"], u_vrnr["u_Xim"], u_vrnr["cov_RI"]),
                (u_vrxr["u_Xre"], u_vrxr["u_Xim"], u_vrxr["cov_RI"])
            ))

        console.print(
            "Factor k expresado en Magnitud y fase.\n"
            "Mag = Rx/Rn\n"
            "Phase = ΦRx - ΦRn "
        )

        k_unc = uncertainty.correction_uncertainty(voltages, uncertainties)

        tables.print_k_table(console, k_values, k_unc, title="Factor k")
    #--------------------------------------------------------------------
    # 10-Calcular Resistor incognito
    #--------------------------------------------------------------------
        console.rule("[bold green]Resistor Incognito")
        progress.update(task, description="[cyan]Calculando Resistor")

        Rx_results = uncertainty.calc_Rx(
        k_values,
        k_unc,
        Rn_mag,
        Rn_pha,
        U_Rn_mag,
        U_Rn_pha
        )
        grapher.add_errorbar(
        x=range(len(Rx_results["Rx_phase"])),
        y=Rx_results["Rx_phase"],
        yerr=Rx_results["u_Rx_phase"],
        titulo="Fase de Rx",
        xlabel="Medición",
        ylabel="Fase [rad]"
    )

    #--------------------------------------------------------------------
    # 11- Promedios
    #--------------------------------------------------------------------
        console.print("|Rx| = |Rn| * |k|\nφ_Rx = φ_Rn + φ_k")     
        console.print(f"f0 = {tools.format_with_prefix(fund_freq, 'Hz')}")
        console.print(f"|Rn| = {tools.format_with_prefix(Rn_mag, 'Ω',5)} ± {tools.format_with_prefix(U_Rn_mag, 'Ω')}")
        console.print(f"φ_Rn = {tools.format_with_prefix(Rn_pha, 'rad',5)} ± {tools.format_with_prefix(U_Rn_pha, 'rad')}")
        tables.print_Rx_table(console, Rx_results, title="Rx",fund_freq=fund_freq)
        
        
    #--------------------------------------------------------------------
        if GRAPH_INFO:
            grapher.show()
        console.rule("[bold green]Hecho")
        progress.update(task, description="[bold green]Listo! ✔")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Procesador de mediciones"
    )

    parser.add_argument(
        "--file",
        type=str,
        default=DEFAULT_MEASURE_PATH
    )

    parser.add_argument(
        "--rn_mag",
        type=float,
        default=DEFAULT_Rn_mag
    )

    parser.add_argument(
        "--u_rn_mag",
        type=float,
        default=DEFAULT_U_Rn_mag
    )

    parser.add_argument(
        "--rn_pha",
        type=float,
        default=DEFAULT_Rn_pha
    )

    parser.add_argument(
        "--u_rn_pha",
        type=float,
        default=DEFAULT_U_Rn_pha
    )

    args = parser.parse_args()

    print(args)

    main(
        Rn_mag=args.rn_mag,
        U_Rn_mag=args.u_rn_mag,
        Rn_pha=args.rn_pha,
        U_Rn_pha=args.u_rn_pha,
        MEASURE_PATH=args.file,
    )