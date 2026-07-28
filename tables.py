from rich.table import Table
import numpy as np
import tools
from rich.panel import Panel

def print_k_table(console, k_values, k_unc, title="Resultados de k"):

    table = Table(title=title)
    table.add_column("Set", justify="right")

    for i in range(len(k_values)):
        table.add_column(f"{i+1}", justify="right")

    table.add_column("Avg", justify="right")
    table.add_column("u Tipo A", justify="right")

    #--------------------------------------
    # Magnitud
    #--------------------------------------
    mags = np.abs(k_values)
    u_mags = k_unc["u_k_mag"]

    mag_cells = [
        f"{tools.format_with_prefix(m, '')} ± {tools.format_with_prefix(u, '')}"
        for m, u in zip(mags, u_mags)
    ]

    table.add_row(
        "Mag",
        *mag_cells,
        tools.format_with_prefix(np.mean(mags), ''),
        tools.format_with_prefix(np.std(mags), '')
    )

    #--------------------------------------
    # Fase
    #--------------------------------------
    phases = np.angle(k_values)
    u_phases = k_unc["u_k_phase"]

    phase_cells = [
        f"{tools.format_with_prefix(p, 'rad')} ± {tools.format_with_prefix(u, 'rad')}"
        for p, u in zip(phases, u_phases)
    ]

    table.add_row(
        "Phase",
        *phase_cells,
        tools.format_with_prefix(np.mean(phases), 'rad'),
        tools.format_with_prefix(np.std(phases), 'rad')
    )

    console.print(table)

def print_fft_table(console, fft_results, fft_uncertainties, title):

    table = Table(title=title)
    table.add_column("Set", justify="right")

    for i in range(len(fft_results)):
        table.add_column(f"{i+1}", justify="right")

    table.add_column("Mean", justify="right")
    table.add_column("Std", justify="right")

    mags = []
    phases = []
    freqs = []

    mags_u = []
    phases_u = []

    for fft, unc in zip(fft_results, fft_uncertainties):

        idx = fft["fundamental_freq_index"]

        X = fft["dft"][idx]
        Re = np.real(X)
        Im = np.imag(X)

        A = np.sqrt(Re**2 + Im**2)
        phi = np.arctan2(Im, Re)
        f0 = fft["freqs"][idx]

        mags.append(A)
        phases.append(phi)
        freqs.append(f0)

        mags_u.append(unc["u_A"])
        phases_u.append(unc["u_phi"])

    # -------------------------
    # Promedios
    # -------------------------

    mean_mag = np.mean(mags)
    std_mag = np.std(mags)

    # Promedio circular de fase
    mean_phase = np.angle(np.mean(np.exp(1j * np.array(phases))))
    std_phase = np.std(phases)  # lo dejás como métrica de dispersión simple

    mean_freq = np.mean(freqs)
    std_freq = np.std(freqs)

    # -------------------------
    # Filas
    # -------------------------

    table.add_row(
        "F0",
        *[tools.format_with_prefix(v, 'Hz') for v in freqs],
        #tools.format_with_prefix(mean_freq, 'Hz'),
        #tools.format_with_prefix(std_freq, 'Hz')
        '--','--'
    )

    table.add_row(
        "Mag",
        *[
            f"{tools.format_with_prefix(v, 'V')} ± {tools.format_with_prefix(u, 'V')}"
            for v, u in zip(mags, mags_u)
        ],
        tools.format_with_prefix(mean_mag, 'V'),
        tools.format_with_prefix(std_mag, 'V')
    )

    table.add_row(
        "Phase",
        *[
            f"{tools.format_with_prefix(v, 'rad')} ± {tools.format_with_prefix(u, 'rad')}"
            for v, u in zip(phases, phases_u)
        ],
        tools.format_with_prefix(mean_phase, 'rad'),
        tools.format_with_prefix(std_phase, 'rad')
    )

    console.print(table)
    
def print_Rx_table(console, Rx_results, title="Resultados de Rx", fund_freq=None):

    table = Table(title=title)
    table.add_column("Set", justify="right")

    n = len(Rx_results["Rx_mag"])

    for i in range(n):
        table.add_column(f"{i+1}", justify="right")

    table.add_column("Avg", justify="right")

    #--------------------------------------
    # Magnitud
    #--------------------------------------
    mags = Rx_results["Rx_mag"]
    u_mags = Rx_results["u_Rx_mag"]

    mag_cells = [
        f"{tools.format_with_prefix(m, 'Ω')} ± {tools.format_with_prefix(u, 'Ω')}"
        for m, u in zip(mags, u_mags)
    ]

    table.add_row(
        "|Rx|",
        *mag_cells,
        f"{tools.format_with_prefix(np.mean(mags), 'Ω')} ± {tools.format_with_prefix(np.mean(u_mags))}",
       
    )

    #--------------------------------------
    # Fase
    #--------------------------------------
    phases = Rx_results["Rx_phase"]
    u_phases = Rx_results["u_Rx_phase"]

    phase_cells = [
        f"{tools.format_with_prefix(p, 'rad')} ± {tools.format_with_prefix(u, 'rad')}"
        for p, u in zip(phases, u_phases)
    ]

    table.add_row(
        "φ_Rx",
        *phase_cells,
        f"{tools.format_with_prefix(np.mean(phases), 'rad')} ± {tools.format_with_prefix(np.mean(u_phases), 'rad')}",
    )

    console.print(table)
    console.print(Panel.fit(f"{tools.format_with_prefix(np.mean(mags), 'Ω')} ± {tools.format_with_prefix(np.mean(u_mags),'Ω')} || {tools.format_with_prefix(np.mean(phases), 'rad')} ± {tools.format_with_prefix(np.mean(u_phases),'rad')} || {tools.format_with_prefix((np.mean(phases)/(np.pi*2*fund_freq)), 's')} ± {tools.format_with_prefix((np.mean(u_phases)/(np.pi*2*fund_freq)),'s')}", style="bold white on green"))
    console.print(
    Panel.fit(
        f"{fund_freq:.1f}".replace(".", ",") + ";" +
        f"{np.mean(mags):.3e}".replace(".", ",") + ";" +
        f"{np.mean(u_mags):.3e}".replace(".", ",") + ";" +
        f"{np.mean(phases):.3e}".replace(".", ",") + ";" +
        f"{np.mean(u_phases):.3e}".replace(".", ",") + ";" +
        f"{(np.mean(phases)/(2*np.pi*fund_freq)):.3e}".replace(".", ",") + ";" +
        f"{(np.mean(u_phases)/(2*np.pi*fund_freq)):.3e}".replace(".", ",")
        #style="bold white on green"
    )
)