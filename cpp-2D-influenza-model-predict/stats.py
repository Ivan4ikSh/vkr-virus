import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.mplot3d import Axes3D  # необходимо для проекции 3d

def parse_args():
    parser = argparse.ArgumentParser(description='Visualize 2D antigenic evolution (separate I and R profiles)')
    parser.add_argument('--data-dir', type=str, required=True, help='Directory with CSV files')
    parser.add_argument('--output-dir', type=str, default='.', help='Directory to save plots')
    parser.add_argument('--max-slices', type=int, default=10, help='Maximum number of time slices to show')
    parser.add_argument('--no-time-series', action='store_true', help='Do not create time series plot')
    return parser.parse_args()

def read_state(file_pattern):
    files = sorted(glob.glob(file_pattern), key=lambda x: int(x.split('_')[-1].split('.')[0]))
    states = []
    steps = []
    for f in files:
        df = pd.read_csv(f, sep=';', header=None)
        states.append(df.values)
        steps.append(int(f.split('_')[-1].split('.')[0]))
    return states, steps

def main():
    args = parse_args()
    data_dir = args.data_dir
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    # --- Загрузка состояний I и R ---
    I_states, I_steps = read_state(os.path.join(data_dir, 'state_I_step_*.csv'))
    R_states, R_steps = read_state(os.path.join(data_dir, 'state_R_step_*.csv'))

    if not I_states:
        print("No state files found. Exiting.")
        return

    # Синхронизация шагов (на случай несовпадения)
    if I_steps != R_steps:
        print("Warning: step mismatch between I and R files. Using intersection.")
        common_steps = sorted(set(I_steps) & set(R_steps))
        I_states = [I_states[I_steps.index(s)] for s in common_steps]
        R_states = [R_states[R_steps.index(s)] for s in common_steps]
        I_steps = common_steps
        R_steps = common_steps

    # Определяем размер сетки
    Ly, Lx = I_states[0].shape  # строки = Y, столбцы = X
    x_axis = np.arange(Lx)
    y_axis = np.arange(Ly)
    X, Y = np.meshgrid(x_axis, y_axis)

    # Выбираем не более max_slices равномерно
    n_slices = len(I_states)
    if n_slices > args.max_slices:
        indices = np.linspace(0, n_slices - 1, args.max_slices, dtype=int)
    else:
        indices = range(n_slices)

    selected_steps = [I_steps[i] for i in indices]
    selected_I = [I_states[i] for i in indices]
    selected_R = [R_states[i] for i in indices]

    # --- Построение 3D поверхностей для I и R ---
    n = len(selected_steps)
    if n == 0:
        print("No slices selected.")
        return

    # Определяем сетку подграфиков (квадратную или почти квадратную)
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))

    # Фигура для I
    fig_I = plt.figure(figsize=(cols * 5, rows * 4))
    for idx, (step, I_mat) in enumerate(zip(selected_steps, selected_I)):
        ax = fig_I.add_subplot(rows, cols, idx + 1, projection='3d')
        surf = ax.plot_surface(X, Y, I_mat, cmap='hot', edgecolor='none', alpha=0.8)
        ax.set_title(f'I, step {step}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('I')
        ax.view_init(elev=30, azim=45)  # фиксированный угол обзора
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'I_3d.png'), dpi=150)
    print(f"3D I plot saved to {os.path.join(out_dir, 'I_3d.png')}")

    # Фигура для R
    fig_R = plt.figure(figsize=(cols * 5, rows * 4))
    for idx, (step, R_mat) in enumerate(zip(selected_steps, selected_R)):
        ax = fig_R.add_subplot(rows, cols, idx + 1, projection='3d')
        surf = ax.plot_surface(X, Y, R_mat, cmap='Blues', edgecolor='none', alpha=0.8)
        ax.set_title(f'R, step {step}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('R')
        ax.view_init(elev=30, azim=45)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'R_3d.png'), dpi=150)
    print(f"3D R plot saved to {os.path.join(out_dir, 'R_3d.png')}")

    # --- Временные ряды (опционально, без изменений) ---
    if not args.no_time_series:
        norm_df = pd.read_csv(os.path.join(data_dir, 'norm_time_series.csv'), sep=';')
        finf_df = pd.read_csv(os.path.join(data_dir, 'finf_time_series.csv'), sep=';')
        steps_ts = norm_df['Step'].values
        norm = norm_df['Norm'].values
        finf = finf_df['Finf'].values

        fig2, ax2 = plt.subplots(figsize=(10, 6))
        #ax2.plot(steps_ts, norm, 'k-', label='Norm (total/N)')
        ax2.plot(steps_ts, finf, 'r-', label='Fraction infected')

        meanx_path = os.path.join(data_dir, 'mean_x_time_series.csv')
        if os.path.exists(meanx_path):
            meanx_df = pd.read_csv(meanx_path, sep=';')
            meanx = meanx_df['MeanX'].values
            ax2.plot(steps_ts, meanx, 'b--', label='Mean X')

        meany_path = os.path.join(data_dir, 'mean_y_time_series.csv')
        if os.path.exists(meany_path):
            meany_df = pd.read_csv(meany_path, sep=';')
            meany = meany_df['MeanY'].values
            ax2.plot(steps_ts, meany, 'g--', label='Mean Y')

        ax2.set_xlabel('Step')
        ax2.set_ylabel('Value')
        ax2.set_title('Time series')
        ax2.legend(loc='best')
        ax2.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'time_series.png'), dpi=150)
        print(f"Time series plot saved to {os.path.join(out_dir, 'time_series.png')}")

    print("Done.")

if __name__ == '__main__':
    main()