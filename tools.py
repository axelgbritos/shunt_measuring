'''
tools.py


V1.0 03/2026

                Axel Britos
Instituto Nacional de Tecnologia Industrial

~~~Notas de version~~~
V1.0 03/2026. Version inicial

'''

def format_with_prefix(value, unit='', sig_digits=4):
    """Formatea un valor utilizando prefijos SI y la cantidad de cifras significativas indicada."""

    prefixes = [
        (1e-9, 'n'),
        (1e-6, 'μ'),
        (1e-3, 'm'),
        (1, ''),
        (1e3, 'k'),
        (1e6, 'M'),
        (1e9, 'G'),
    ]

    abs_val = abs(value)

    for factor, prefix in reversed(prefixes):
        if abs_val >= factor:
            return f"{value/factor:.{sig_digits}g}{prefix}{unit}"

    return f"{value:.{sig_digits}g}{unit}"


#TODO: Funciones para la lectura de variables desde la ejecucion del script con argumentos incluidos.
'''
parser = argparse.ArgumentParser(description="Procesamiento de mediciones")

parser.add_argument(
    "file",
    nargs="?",                 # hace que el argumento sea opcional
    default=DEFAULT_PATH,
    help="archivo de medición (csv)"
)

args = parser.parse_args()
'''