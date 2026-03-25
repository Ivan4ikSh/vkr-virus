import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob
import re
import argparse
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

class VirusWaveVisualizer2D:
    def __init__(self, data_dir="out/data"):
        """
        Visualizer for the 2D model.
        data_dir : folder containing I_step_*.csv, R_step_*.csv, finf.csv and parameters.txt
        """
        self.data_dir = data_dir
        if not os.path.exists(data_dir):
            print(f"Error: folder {data_dir} does not exist!")
            return

        self.find_data_files()
        self.model_params = self.load_model_parameters()

        # Colormap for composite image ("comet")
        colors = ['#FFFFFF', '#87CEEB', '#32CD32', '#FFD700', '#FF4500', '#8B0000']
        self.wave_cmap = LinearSegmentedColormap.from_list('wave', colors, N=256)

    def find_data_files(self):
        """Find all I_step_*.csv and R_step_*.csv files, sort by step number."""
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

    def load_model_parameters(self):
        param_file = os.path.join(self.data_dir, "parameters.txt")
        params = {}
        if os.path.exists(param_file):
            with open(param_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or '=' not in line: continue
                    key, val = line.split('=', 1)
                    key, val = key.strip(), val.strip()
                    try:
                        params[key] = float(val) if ('.' in val or 'e' in val.lower()) else int(val)
                    except ValueError:
                        params[key] = val
            print(f"Loaded {len(params)} parameters")
        else:
            print(f"Parameter file not found, using defaults")
            params = {'L': 0, 'M': 0, 'R0': 0, 'Ub': 0, 'a': 0, 'dt': 0, 'N': 0}
        return params

    def _format_parameters(self):
        display_names = {
            'L': 'Grid size L', 'M': 'Total steps M', 'R0': 'Reproduction number R₀',
            'Ub': 'Mutation rate Ub', 'a': 'Immunity scale a', 'dt': 'Time step dt',
            'init_infected': 'Initial infected', 'N': 'Population N', 'Seed': 'Seed'
        }
        order = ['L', 'M', 'R0', 'Ub', 'a', 'dt', 'init_infected', 'N', 'Seed']
        lines = []
        for key in order:
            if key in self.model_params:
                val = self.model_params[key]
                val_str = f"{val:.2e}" if isinstance(val, float) and (abs(val) < 1e-4 or abs(val) > 1e4) else str(val)
                lines.append(f"{display_names.get(key, key)}: {val_str}")
        return "\n".join(lines)

    def create_finf_plot(self, save_path="finf_evolution.png"):
        """
        Reads finf.csv and creates a diagnostic plot matching the provided example style.
        """
        csv_path = os.path.join(self.data_dir, "finf.csv")
        if not os.path.exists(csv_path):
            print(f"Warning: {csv_path} not found. Skipping finf plot.")
            return

        # 1. Load data
        df = pd.read_csv(csv_path, sep=';')
    
        # 2. Prepare parameter string for title (as in example)
        p = self.model_params
        # We take keys that are present in your C++ struct and display them
        param_str = f"R0={p.get('R0', 'N/A')}, Ub={p.get('Ub', 'N/A')}, a={p.get('a', 'N/A')}, " \
                    f"L={p.get('L', 'N/A')}, N={p.get('N', 'N/A'):.0e}, dt={p.get('dt', 'N/A')}"
    
        # Calculation of some observables (optional, based on your screenshot)
        # If you don't calculate speed/sdl here, we can just show the mean finf
        avg_finf = df['finf'].mean()
        obs_str = f"mean_finf={avg_finf:.2e}"

        # 3. Plotting
        fig, ax = plt.subplots(figsize=(12, 7)) # Adjust size to match aspect ratio
    
        ax.plot(df['step'], df['finf'], 'r-', lw=1.5) # 'r-' is red solid line
    
        # 4. Styling (Matching the example)
        ax.set_xlabel('Time', fontsize=11)
        ax.set_ylabel('Fraction infected (finf)', fontsize=11)
    
        # Double-layered title
        full_title = f"Fraction infected over time\nParams: {param_str}\nObservables: {obs_str}"
        ax.set_title(full_title, fontsize=12, pad=15)
    
        # Dotted grid as in 1D_diagnostics.png
        ax.grid(True, linestyle=':', alpha=0.6)
    
        # Setting limits to look cleaner
        ax.set_xlim(df['step'].min() - 50, df['step'].max() + 50)

        plt.tight_layout()
    
        # 5. Save
        plt.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"Diagnostic plot saved as {save_path}")

    def create_wave_snapshot(self, step_indices=None, save_path="wave_snapshots.png"):
        """
        Create a montage of 4 time steps.
        Top row – infected (I), bottom – susceptible (R).
        Right side – parameter panel.
        """
        if step_indices is None:
            total = len(self.I_files)
            if total == 0:
                print("No data for snapshot")
                return
            # Choose 4 evenly spaced steps
            step_indices = [0, total//3, 2*total//3, total-1]
            step_indices = [max(0, min(idx, total-1)) for idx in step_indices]

        # Reduced figure height: from (20,10) to (16,6)
        fig = plt.figure(figsize=(20, 6))
        gs = gridspec.GridSpec(2, 5, width_ratios=[1,1,1,1,0.8], wspace=0.25, hspace=0.3)

        for col, idx in enumerate(step_indices[:4]):
            if idx >= len(self.I_files):
                continue

            I = self.load_matrix(self.I_files[idx])
            R = self.load_matrix(self.R_files[idx])
            step_num = self.extract_step_number(self.I_files[idx])

            # Infected (top)
            ax_I = plt.subplot(gs[0, col])
            im_I = ax_I.imshow(I, cmap='Reds', aspect='auto', origin='lower', vmin=0)
            ax_I.set_title(f'Step {step_num}', fontsize=10, fontweight='bold')
            if col == 0:
                ax_I.set_ylabel('y (antigenic)', fontsize=8)
            ax_I.tick_params(labelsize=7)
            plt.colorbar(im_I, ax=ax_I, fraction=0.046, pad=0.04)

            # Susceptible (bottom)
            ax_R = plt.subplot(gs[1, col])
            im_R = ax_R.imshow(R, cmap='Greens', aspect='auto', origin='lower', vmin=0)
            ax_R.set_xlabel('x (antigenic)', fontsize=8)
            if col == 0:
                ax_R.set_ylabel('y (antigenic)', fontsize=8)
            ax_R.tick_params(labelsize=7)
            plt.colorbar(im_R, ax=ax_R, fraction=0.046, pad=0.04)

        # Parameter panel on the right (spans both rows)
        ax_params = plt.subplot(gs[:, 4])
        ax_params.axis('off')
        param_text = self._format_parameters()
        ax_params.text(0.05, 0.5, param_text, transform=ax_params.transAxes, fontsize=8, fontfamily='monospace', verticalalignment='center', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

        plt.suptitle('Wave evolution', fontsize=12, fontweight='bold')
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Snapshot saved as {save_path}")

    def create_wave_animation(self, save_path="wave_evolution.gif", max_frames=300):
        """
        Creates an animation of wave evolution with 3 plots and parameters.
        (Unchanged from original)
        """
        if len(self.I_files) == 0:
            print("No data for animation!")
            return
            
        # Limit number of frames
        total_frames = min(max_frames, len(self.I_files))
        
        fig = plt.figure(figsize=(18, 10))
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
        
        # Load first frame for initialization
        I0 = self.load_matrix(self.I_files[0])
        R0 = self.load_matrix(self.R_files[0])
        
        # 1. Composite image (top left)
        ax1 = plt.subplot(gs[0, 0])
        composite0 = I0 * 0.8 + R0 * 0.4
        im1 = ax1.imshow(composite0, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian', animated=True)
        ax1.set_title('Wave evolution (comet)', fontsize=12, fontweight='bold')
        ax1.set_xlabel('Antigenic coordinate (x)')
        ax1.set_ylabel('Antigenic coordinate (y)')
        ax1.grid(True, alpha=0.3, linestyle='--')
        plt.colorbar(im1, ax=ax1, label='Intensity')
        
        # 2. Infected (top right)
        ax2 = plt.subplot(gs[0, 1])
        im2 = ax2.imshow(I0, cmap='Reds', aspect='auto', origin='lower', animated=True)
        ax2.set_title('Infected', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Antigenic coordinate (x)')
        ax2.set_ylabel('Antigenic coordinate (y)')
        plt.colorbar(im2, ax=ax2, label='Infected density')
        
        # 3. Susceptible (bottom left)
        ax3 = plt.subplot(gs[1, 1])
        im3 = ax3.imshow(R0, cmap='Greens', aspect='auto', origin='lower', animated=True)
        ax3.set_title('Recovered', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Antigenic coordinate (x)')
        ax3.set_ylabel('Antigenic coordinate (y)')
        plt.colorbar(im3, ax=ax3, label='Susceptible density')
        
        # 4. Model parameters (bottom right)
        ax4 = plt.subplot(gs[1, 0])
        ax4.axis('off')
        params_text = self._format_parameters()
        ax4.text(0.05, 0.5, params_text, transform=ax4.transAxes, fontsize=15, fontfamily='monospace', verticalalignment='center', horizontalalignment='left', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='gray'))
        
        step_num = self.extract_step_number(self.I_files[0])
        images = [im1, im2, im3]
        
        def update(frame):
            # Load data for current frame
            I = self.load_matrix(self.I_files[frame])
            R = self.load_matrix(self.R_files[frame])
            
            # Update images
            composite = I * 10 + R
            images[0].set_data(composite)
            images[1].set_data(I)
            images[2].set_data(R)
            
            # Autoscale
            for im in images:
                im.autoscale()
            
            # Update step info
            step_num = self.extract_step_number(self.I_files[frame])            
            # Update suptitle
            fig.suptitle(f'Virus wave evolution animation - Step: {step_num}', fontsize=18, fontweight='bold')
            
            return images
        
        # Create animation
        anim = FuncAnimation(fig, update, frames=total_frames, interval=300, blit=True, repeat=True)
        
        print(f"Creating animation with {total_frames} frames...")
        
        # Save animation
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
    if not viz.I_files: return

    viz.create_wave_snapshot(save_path=args.snapshot)
    viz.create_finf_plot(save_path=args.finf_plot)
    viz.create_wave_animation(save_path=args.animation)
    print("DONE")

if __name__ == "__main__":
    main()