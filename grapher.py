import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from matplotlib.gridspec import GridSpec
import numpy as np
import tools

class pageGraphs:
    def __init__(self):
        self._paginas = []  # lista de dicts: {"titulo", "fn", "tipo"}

    # ── Métodos para registrar gráficos ──────────────────────────────────

    def signal_segments(self, signal, segment_size, transient_length, signal_length, titulo="Signal segmentada"):
        def _fn(fig, ax):
            ax.plot(signal, color='steelblue', linewidth=0.8)
            
            n_segments = len(signal) // segment_size
            
            for i in range(n_segments):
                inicio = i * segment_size
                fin_transitorio = inicio + transient_length
                fin_valido = fin_transitorio + signal_length
                
                # Línea roja: inicio de segmento
                ax.axvline(x=inicio, color='red', linewidth=1, linestyle='--')
                # Línea azul: fin del transitorio (inicio de zona válida)
                ax.axvline(x=fin_transitorio, color='blue', linewidth=1, linestyle='--')
                # Línea azul: fin de zona válida
                if fin_valido <= len(signal):
                    ax.axvline(x=fin_valido, color='blue', linewidth=1, linestyle='--')
            
            ax.set_title(titulo)
            ax.set_xlabel("Muestra")
            ax.set_ylabel("Amplitud")
            ax.legend(["Señal de entrada", "Segmento", "Segmento sin transitorio"], loc="upper right")
    
        self._paginas.append({"titulo": titulo, "fn": _fn, "tipo": "simple"})

    def add(self, fn, titulo=""):
        """Registra una función genérica fn(fig, ax) o fn(fig) que dibuja una página."""
        self._paginas.append({"titulo": titulo, "fn": fn, "tipo": "simple"})

    def add_hist(self, datos, titulo="", xlabel="", ylabel="Frecuencia", bins=50):
        def _fn(fig, ax):
            ax.hist(datos, bins=bins)
            ax.set_title(titulo)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        self._paginas.append({"titulo": titulo, "fn": _fn, "tipo": "simple"})

    def add_dispersion(self, x, y, titulo="", xlabel="", ylabel="", **scatter_kwargs):
        def _fn(fig, ax):
            ax.scatter(x, y, **scatter_kwargs)
            ax.set_title(titulo)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        self._paginas.append({"titulo": titulo, "fn": _fn, "tipo": "simple"})
    
    def add_errorbar(self, x, y, yerr=None, xerr=None, titulo="", xlabel="", ylabel="", **errorbar_kwargs):
        def _fn(fig, ax):
            ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o', capsize=5, **errorbar_kwargs)
            ax.axhline(y=np.mean(y), color='red', linewidth=1, linestyle='--', label=f"fase: {tools.format_with_prefix(np.mean(y), unit='rad')} ± {tools.format_with_prefix(np.mean(yerr), unit='rad')}")
            ax.legend()
            ax.set_title(titulo)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        self._paginas.append({"titulo": titulo, "fn": _fn, "tipo": "simple"})

    def add_subplots(self, fns_matriz, titulo="", compartir_x=False, compartir_y=False):
        """
        fns_matriz: lista de listas de funciones fn(ax), define la grilla.
        Ejemplo 2x2: [[fn1, fn2], [fn3, fn4]]
        """
        filas = len(fns_matriz)
        cols  = max(len(fila) for fila in fns_matriz)

        def _fn(fig, _ax_ignorado):
            # Borramos el ax placeholder y creamos la grilla real
            _ax_ignorado.set_visible(False)
            gs = GridSpec(filas, cols, figure=fig,
                          top=0.88, bottom=0.18, hspace=0.45, wspace=0.35)
            for r, fila in enumerate(fns_matriz):
                for c, fn_subax in enumerate(fila):
                    ax_sub = fig.add_subplot(gs[r, c])
                    fn_subax(ax_sub)
            if titulo:
                fig.suptitle(titulo, fontsize=13, y=0.97)

        self._paginas.append({"titulo": titulo, "fn": _fn, "tipo": "subplots"})

    # ── Mostrar todo ──────────────────────────────────────────────────────

    def show(self):
        if not self._paginas:
            print("No hay páginas registradas.")
            return

        fig, ax_main = plt.subplots()
        fig.subplots_adjust(bottom=0.18, top=0.92)
        pagina_actual = [0]

        def dibujar(i):
            # Limpia axes extra (de subplots previos) y el ax principal
            for a in fig.axes:
                if a not in (ax_prev_btn_ax, ax_next_btn_ax):
                    a.cla()
            # Elimina axes adicionales dejando solo el principal
            for a in fig.axes[:]:
                if a not in (ax_main, ax_prev_btn_ax, ax_next_btn_ax):
                    a.remove()

            ax_main.set_visible(True)
            pagina = self._paginas[i]
            pagina["fn"](fig, ax_main)

            # Contador de página
            fig.texts = [t for t in fig.texts if not getattr(t, "_es_contador", False)]
            txt = fig.text(0.5, 0.01, f"{i+1} / {len(self._paginas)}",
                           ha="center", va="bottom", fontsize=9, color="gray")
            txt._es_contador = True
            fig.canvas.draw_idle()

        def next_graph(event):
            pagina_actual[0] = (pagina_actual[0] + 1) % len(self._paginas)
            dibujar(pagina_actual[0])

        def prev_graph(event):
            pagina_actual[0] = (pagina_actual[0] - 1) % len(self._paginas)
            dibujar(pagina_actual[0])

        ax_prev_btn_ax = plt.axes([0.30, 0.05, 0.15, 0.07])
        ax_next_btn_ax = plt.axes([0.55, 0.05, 0.15, 0.07])
        btn_prev = widgets.Button(ax_prev_btn_ax, "← Anterior")
        btn_next = widgets.Button(ax_next_btn_ax, "Siguiente →")
        btn_prev.on_clicked(prev_graph)
        btn_next.on_clicked(next_graph)

        dibujar(0)
        plt.show()