import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob
import re
import argparse
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

BASE_SIZE = 14
plt.rcParams.update({
    'font.size': BASE_SIZE,
    'axes.titlesize': BASE_SIZE + 2,
    'axes.labelsize': BASE_SIZE + 2,
    'xtick.labelsize': BASE_SIZE - 2,
    'ytick.labelsize': BASE_SIZE - 2,
    'legend.fontsize': BASE_SIZE - 2,
    'figure.titlesize': BASE_SIZE + 2,
    'figure.dpi': 150
})

class VirusWaveVisualizer2D:
    def __init__(self, data_dir="out/data"):
        self.data_dir = data_dir
        if not os.path.exists(data_dir):
            print(f"Error: folder {data_dir} does not exist!")
            return

        self.find_data_files()
        self.params_title_str = self.load_model_parameters_string()

        # --- РАСЧЕТ КОЭФФИЦИЕНТА ДЛЯ ПЕРЕВОДА В ГОДЫ ---
        DAYS_IN_CYCLE = 5.2
        DAYS_IN_YEAR = 365.25
        self.dt = 0.5  # Значение по умолчанию
        
        # Пытаемся вытащить dt из строки параметров с помощью регулярного выражения
        match = re.search(r'dt\s*=\s*([0-9.]+)', self.params_title_str)
        if match:
            self.dt = float(match.group(1))
            
        # Формула: (шаги * dt) = циклы -> (циклы * 5) / 365.25 = годы
        self.TO_YEARS = (self.dt * DAYS_IN_CYCLE) / DAYS_IN_YEAR
        # -----------------------------------------------

        colors = ['#FFFFFF', '#87CEEB', '#32CD32', '#FFD700', '#FF4500', '#8B0000']
        self.wave_cmap = LinearSegmentedColormap.from_list('wave', colors, N=256)

    def find_data_files(self):
        pattern_I = os.path.join(self.data_dir, "I_step_*.csv")
        pattern_R = os.path.join(self.data_dir, "R_step_*.csv")
        self.I_files = sorted(glob.glob(pattern_I), key=self.extract_step_number)
        self.R_files = sorted(glob.glob(pattern_R), key=self.extract_step_number)

        if self.I_files:
            print(f"Found {len(self.I_files)} steps in {self.data_dir}")
        else:
            print(f"No I_step_*.csv files found in {self.data_dir}")

    @staticmethod
    def extract_step_number(filename):
        match = re.search(r'_step_(\d+)', filename)
        return int(match.group(1)) if match else 0

    @staticmethod
    def load_matrix(filename):
        return np.loadtxt(filename, delimiter=';')

    def load_model_parameters_string(self):
        param_file = os.path.join(self.data_dir, "parameters.txt")
        if os.path.exists(param_file):
            with open(param_file, 'r') as f:
                return f.read().strip()
        else:
            print("Parameter file not found!")
            return "Parameters: N/A"

    def create_finf_plot(self, save_path="finf_evolution.png"):
        csv_path = os.path.join(self.data_dir, "finf.csv")
        if not os.path.exists(csv_path):
            print(f"Warning: {csv_path} not found. Skipping finf plot.")
            return

        df = pd.read_csv(csv_path, sep=';')
        avg_finf = df['finf'].mean()
        obs_str = f"mean_finf={avg_finf:.2e}"

        # --- ПРИМЕНЯЕМ ФОРМУЛУ ПЕРЕВОДА В ГОДЫ ---
        df['year'] = df['step'] * self.TO_YEARS

        fig, ax = plt.subplots(figsize=(14, 9))
        
        # Строим график от новой колонки 'year'
        ax.plot(df['year'], df['finf'], 'r-', lw=2)
    
        ax.set_xlabel('Time (years)') # Изменили подпись
        ax.set_ylabel('Fraction infected (finf)')
    
        full_title = f"Fraction infected over time\n\n{self.params_title_str}\n\nObservables: {obs_str}"
        ax.set_title(full_title, pad=20)
    
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Обновляем пределы графика по оси X
        ax.set_xlim(df['year'].min(), df['year'].max())

        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig)
        print(f"Diagnostic plot saved as {save_path}")

    def create_wave_snapshot(self, step_indices=None, save_path="wave_snapshots.png"):
        if step_indices is None:
            total = len(self.I_files)
            if total == 0: return
            step_indices = [0, total//3, 2*total//3, total-1]
            step_indices = [max(0, min(idx, total-1)) for idx in step_indices]

        fig = plt.figure(figsize=(24, 11))
        gs = gridspec.GridSpec(2, 4, wspace=0.075, hspace=0.075)

        for col, idx in enumerate(step_indices[:4]):
            I = self.load_matrix(self.I_files[idx])
            R = self.load_matrix(self.R_files[idx])
            step_num = self.extract_step_number(self.I_files[idx])
            
            # Считаем текущий год для заголовка
            year_val = step_num * self.TO_YEARS

            # --- Ряд Infected ---
            ax_I = plt.subplot(gs[0, col])
            im_I = ax_I.imshow(I, cmap='Reds', aspect='auto', origin='lower', vmin=0)
            ax_I.set_title(f'Year {year_val:.1f}', fontweight='bold') # Подписываем годами
            ax_I.set_xticks([])
            if col == 0:
                ax_I.set_ylabel('Infected (I)\ny (neutral)')
            
            cbar_I = plt.colorbar(im_I, ax=ax_I)
            cbar_I.ax.tick_params(labelsize=BASE_SIZE)
        
            fmt_I = ticker.ScalarFormatter(useMathText=False)
            fmt_I.set_scientific(True)
            fmt_I.set_powerlimits((0, 0))
            cbar_I.ax.yaxis.set_major_formatter(fmt_I)

            # --- Ряд Recovered ---
            ax_R = plt.subplot(gs[1, col])
            im_R = ax_R.imshow(R, cmap='Greens', aspect='auto', origin='lower', vmin=0)
            ax_R.set_xlabel('x (antigenic)')
            if col == 0:
                ax_R.set_ylabel('Recovered (R)\ny (neutral)')
                
            cbar_R = plt.colorbar(im_R, ax=ax_R)
            cbar_R.ax.tick_params(labelsize=BASE_SIZE)
        
            fmt_R = ticker.ScalarFormatter(useMathText=False)
            fmt_R.set_scientific(True)
            fmt_R.set_powerlimits((0, 0))
            cbar_R.ax.yaxis.set_major_formatter(fmt_R)

        fig.suptitle(f'Wave Evolution\n{self.params_title_str}', fontweight='bold', y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 1])

        plt.savefig(save_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Snapshot saved as {save_path}")

    def create_wave_animation(self, save_path="wave_evolution.gif", max_frames=300):
        if len(self.I_files) == 0:
            print("No data for animation!")
            return
            
        total_frames = min(max_frames, len(self.I_files))
        
        fig = plt.figure(figsize=(22, 9))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
        
        I0 = self.load_matrix(self.I_files[0])
        R0 = self.load_matrix(self.R_files[0])
        
        # Composite image
        ax1 = plt.subplot(gs[0])
        composite0 = I0 * 0.8 + R0 * 0.4
        im1 = ax1.imshow(composite0, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian', animated=True)
        ax1.set_title('Wave evolution (comet)', fontweight='bold')
        ax1.set_xlabel('Antigenic (x)')
        ax1.set_ylabel('Antigenic (y)')
        ax1.grid(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        plt.colorbar(im1, ax=ax1)
        
        # Infected
        ax2 = plt.subplot(gs[1])
        im2 = ax2.imshow(I0, cmap='Reds', aspect='auto', origin='lower', animated=True)
        ax2.set_title('Infected', fontweight='bold')
        ax2.set_xlabel('Antigenic (x)')
        ax2.set_ylabel('Antigenic (y)')
        plt.colorbar(im2, ax=ax2)
        
        # Susceptible
        ax3 = plt.subplot(gs[2])
        im3 = ax3.imshow(R0, cmap='Greens', aspect='auto', origin='lower', animated=True)
        ax3.set_title('Recovered', fontweight='bold')
        ax3.set_xlabel('Antigenic (x)')
        ax3.set_ylabel('Antigenic (y)')
        plt.colorbar(im3, ax=ax3)
        
        images = [im1, im2, im3]

        step_num_0 = self.extract_step_number(self.I_files[0])
        year_0 = step_num_0 * self.TO_YEARS
        fig.suptitle(f'Virus wave evolution - Year: {year_0:.1f}\n{self.params_title_str}', fontweight='bold')
        
        plt.tight_layout(rect=[0, 0, 1, 0.90])
        
        def update(frame):
            I = self.load_matrix(self.I_files[frame])
            R = self.load_matrix(self.R_files[frame])
            
            composite = I * 10 + R
            images[0].set_data(composite)
            images[1].set_data(I)
            images[2].set_data(R)
            
            for im in images:
                im.autoscale()
            
            step_num = self.extract_step_number(self.I_files[frame]) 
            year_val = step_num * self.TO_YEARS
            fig.suptitle(f'Virus wave evolution - Year: {year_val:.1f}\n{self.params_title_str}', fontweight='bold')
            
            return images
        
        anim = FuncAnimation(fig, update, frames=total_frames, interval=300, blit=False, repeat=True)
        print(f"Creating animation with {total_frames} frames...")
        
        if save_path.endswith('.gif'):
            anim.save(save_path, writer=PillowWriter(fps=10))
        else:
            anim.save(save_path, writer='ffmpeg', fps=10)
        
        print(f"Animation saved as {save_path}")
        return anim

def main():
    parser = argparse.ArgumentParser(description='2D Epidemic Simulation Visualizer')
    parser.add_argument('-d', '--data-dir', default='out/data')
    parser.add_argument('-a', '--animation', default='wave_evolution.gif')
    parser.add_argument('-s', '--snapshot', default='wave_snapshots.png')
    parser.add_argument('-f', '--finf-plot', default='infection_dynamics.png')
    args = parser.parse_args()

    viz = VirusWaveVisualizer2D(args.data_dir)
    if not hasattr(viz, 'I_files') or not viz.I_files: 
        return

    viz.create_wave_snapshot(save_path=args.snapshot)
    viz.create_finf_plot(save_path=args.finf_plot)
    viz.create_wave_animation(save_path=args.animation)
    print("DONE")

if __name__ == "__main__":
    main()
