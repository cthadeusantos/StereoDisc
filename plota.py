import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def save_heatmap(matrix, filepath, title="Heatmap", cmap="viridis",
                 figsize=(8, 6), show_values=False, fmt=".2f",
                 colorbar=True, dpi=150, show=False, x_tick_labels=None):
    """
    Plota um heatmap de uma matriz (valores entre 0 e 1) e salva em arquivo.
    
    Parâmetros:
    -----------
    matrix : np.ndarray
        Matriz quadrada com valores entre 0 e 1.
    filepath : str ou Path
        Caminho do arquivo de saída (ex: 'heatmap.png', 'figuras/hm.pdf').
    title : str
        Título do gráfico.
    cmap : str
        Colormap do matplotlib (ex: 'viridis', 'hot', 'coolwarm', 'Blues').
    figsize : tuple
        Tamanho da figura em polegadas.
    show_values : bool
        Se True, exibe os valores numéricos em cada célula.
    fmt : str
        Formato dos valores exibidos (usado se show_values=True).
    colorbar : bool
        Se True, exibe a barra de cores.
    dpi : int
        Resolução do arquivo salvo.
    show : bool
        Se True, também exibe a figura na tela (útil em notebooks).
    
    Retorna:
    --------
    Path : caminho do arquivo salvo.
    """
    # Validações
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError("A matriz deve ser bidimensional.")
    finite_values = matrix[np.isfinite(matrix)]
    if finite_values.size and not np.all((finite_values >= 0) & (finite_values <= 1)):
        raise ValueError("Todos os valores devem estar entre 0 e 1.")
    
    n = matrix.shape[0]
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    masked_matrix = np.ma.masked_invalid(matrix)
    display_cmap = plt.get_cmap(cmap).copy()
    display_cmap.set_bad(color="#bdbdbd")
    im = ax.imshow(masked_matrix, cmap=display_cmap, vmin=0, vmax=1,
                   interpolation="nearest")
    has_undefined = np.any(~np.isfinite(matrix))
    if matrix.shape[0] == 1:
        # Stretch the single comparison row to use the compact figure height.
        ax.set_aspect("auto")
    
    if x_tick_labels is None:
        # Sparse ticks keep labels readable for datasets with many trees.
        tick_count = min(matrix.shape[1], 12)
        tick_positions = np.unique(np.linspace(
            0, matrix.shape[1] - 1, tick_count, dtype=int
        )) if matrix.shape[1] else np.array([], dtype=int)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_positions + 1)
    else:
        if len(x_tick_labels) != 2:
            raise ValueError("x_tick_labels deve conter exatamente dois rótulos.")
        if matrix.shape[1] == 1:
            ax.set_xticks([0])
            label = x_tick_labels[0] if x_tick_labels[0] == x_tick_labels[1] else f"{x_tick_labels[0]} / {x_tick_labels[1]}"
            ax.set_xticklabels([label])
        else:
            ax.set_xticks([0, matrix.shape[1] - 1])
            ax.set_xticklabels(x_tick_labels)
    if n == 1:
        ax.set_yticks([])
        ax.set_ylabel("")
    else:
        y_tick_count = min(n, 12)
        y_tick_positions = np.unique(np.linspace(0, n - 1, y_tick_count, dtype=int))
        ax.set_yticks(y_tick_positions)
        ax.set_yticklabels(y_tick_positions + 1)
        ax.set_ylabel("Metric")
    ax.set_xlabel("Tree" if x_tick_labels is None else "")
    ax.set_title(title)
    if has_undefined:
        ax.text(0.01, 1.02, "Cinza: sem pares comparáveis", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=8, color="#555555")
    
    if show_values:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix[i, j]
                if not np.isfinite(val):
                    continue
                color = "white" if val < 0.5 else "black"
                ax.text(j, i, format(val, fmt),
                        ha="center", va="center", color=color, fontsize=9)
    
    if colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.025)
        cbar.set_label("Value")
    
    fig.tight_layout()
    
    # SALVA o arquivo (bbox_inches='tight' evita cortes nas bordas)
    fig.savefig(filepath, dpi=dpi, bbox_inches="tight")
    
    if show:
        plt.show()
    
    # Libera a memória da figura (importante em loops!)
    plt.close(fig)
    
    return filepath



def get_filename_with_new_extension(filepath: str, new_extension: str, mode: str =None) -> str:
    """
    Recebe um caminho de arquivo .nwk ou .newick (com ou sem diretório)
    e retorna apenas o nome do arquivo (sem diretório) com uma nova extensão.
    
    Exemplos:
        'nome_do_arquivo.nwk'                       -> 'nome_do_arquivo.csv'
        './nome_do_arquivo.nwk'                     -> 'nome_do_arquivo.csv'
        'datasets/sintetics/nome_do_arquivo.newick' -> 'nome_do_arquivo.csv'
        'foo/bar/tree.nwk' + '.txt'                 -> 'tree.txt'
    
    Parâmetros:
    -----------
    filepath : str
        Caminho do arquivo de entrada (contém a extensão .nwk ou .newick).
    new_extension : str
        Nova extensão desejada, com ou sem ponto (ex: 'csv', '.txt', 'json').
    
    Retorna:
    --------
    str : nome do arquivo sem diretório e com a nova extensão.
    """
    path = Path(filepath)
    
    # Valida se a extensão original é suportada
    if path.suffix.lower() not in {".nwk", ".newick"}:
        raise ValueError(
            f"Extensão não suportada: '{path.suffix}'. "
            f"Esperado '.nwk' ou '.newick'."
        )
    
    # Normaliza a nova extensão (garante que comece com '.')
    new_ext = new_extension if new_extension.startswith(".") else f".{new_extension}"
    
    # path.stem = nome sem diretório e sem a extensão antiga
    if mode is None:
        mode = "unknown"
    return f"{path.stem}_heatmap_{mode}{new_ext}"

# ----------------- Exemplo de uso -----------------
if __name__ == "__main__":
    np.random.seed(42)
    m = np.random.rand(5, 5)
    
    # Salva como PNG
    save_heatmap(m, "saida/heatmap.png", title="Matriz Aleatória 5x5",
                 cmap="viridis", show_values=True)
    
    # Salva como PDF (vetorial, ótimo para relatórios/artigos)
    save_heatmap(m, "saida/heatmap.pdf", title="Matriz Aleatória 5x5",
                 cmap="hot", dpi=300)
