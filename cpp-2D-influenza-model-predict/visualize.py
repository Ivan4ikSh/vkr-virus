import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob
import re
import argparse
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.gridspec as gridspec

class VirusWaveVisualizer:
    def __init__(self, data_dir="output"):
        """
        Initializes the visualizer for virus evolution wave model.
        Args:
            data_dir (str): Directory containing data files
        """
        self.data_dir = data_dir
        
        if not os.path.exists(data_dir):
            print(f"Error: Folder {data_dir} does not exist!")
            return
            
        self.find_data_files()
        
        # Create custom colormap for "comet" display
        colors = ['#FFFFFF', '#87CEEB', '#32CD32', '#FFD700', '#FF4500', '#8B0000']
        self.wave_cmap = LinearSegmentedColormap.from_list('wave', colors, N=256)
        
        # Load model parameters from file
        self.model_params = self.load_model_parameters()
        
    def load_model_parameters(self):
        """Loads model parameters from parameters.txt file."""
        param_file = os.path.join(self.data_dir, "parameters.txt")
        
        if not os.path.exists(param_file):
            print(f"Parameter file not found: {param_file}")
            # Return default parameters
            return self.get_default_parameters()
        
        try:
            params = {}
            current_section = None
            
            with open(param_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('==='):
                        continue
                    
                    if line.endswith(':'):
                        current_section = line[:-1]
                        continue
                    
                    # Parse lines like "Parameter: value"
                    if ': ' in line:
                        key, value = line.split(': ', 1)
                        
                        # Convert numeric values
                        if key in ['Grid size L', 'Reproductive number R0', 
                                  'Mutation probability D', 'Immunity distance a',
                                  'Coordinate spread Varx', 'Initial infected fraction',
                                  'Show interval tshow', 'Start time T0', 'Total steps M']:
                            try:
                                # Handle scientific notation
                                if 'e' in value.lower() or 'E' in value.lower():
                                    params[key] = float(value)
                                else:
                                    # Check if integer
                                    if '.' not in value and value.replace('.', '').isdigit():
                                        params[key] = int(value)
                                    else:
                                        params[key] = float(value)
                            except ValueError:
                                params[key] = value
                        else:
                            params[key] = value
            
            # Add additional information
            if len(self.I_files) > 0:
                params['Total simulation files'] = len(self.I_files)
                if len(self.I_files) > 0:
                    params['First step'] = self.extract_step_number(self.I_files[0])
                    params['Last step'] = self.extract_step_number(self.I_files[-1])
                    
                    # Load one file to determine grid size
                    sample_data = self.load_matrix(self.I_files[0])
                    params['Data dimensions'] = f"{sample_data.shape[0]}x{sample_data.shape[1]}"
            
            print(f"✓ Loaded {len(params)} model parameters")
            return params
            
        except Exception as e:
            print(f"Error loading parameters: {e}")
            return self.get_default_parameters()
    
    def get_default_parameters(self):
        """Returns default parameters if file is not found."""
        return {
            "Grid size L": "50x50",
            "Reproductive number R0": 2.6,
            "Mutation probability D": 0.0001,
            "Immunity distance a": 7.0,
            "Coordinate spread Varx": 0.01,
            "Total population N": 1.0e10,
            "Initial infected fraction": 0.01,
            "Show interval tshow": 100,
            "Start time T0": 0,
            "Total steps M": 1500,
            "Asymmetry type": 1.0,
            "Note": "Parameters loaded from default values"
        }
    
    def find_data_files(self):
        """Finds all result files in the output folder."""
        self.I_files = sorted(glob.glob(os.path.join(self.data_dir, "state_I_step_*.csv")), key=self.extract_step_number)
        self.R_files = sorted(glob.glob(os.path.join(self.data_dir, "state_R_step_*.csv")), key=self.extract_step_number)
        
        if len(self.I_files) > 0:
            print(f"Found {len(self.I_files)} simulation steps in directory '{self.data_dir}'")
        else:
            print(f"No data files found in directory '{self.data_dir}'")
        
    def extract_step_number(self, filename):
        """Extracts step number from filename."""
        match = re.search(r'step_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    def load_matrix(self, filename):
        """Loads matrix from CSV file."""
        return np.loadtxt(filename, delimiter=';')

    def create_wave_snapshot(self, step_indices=None, save_path="wave_snapshots.png"):
        """
        Creates snapshots of the wave at different stages of evolution.
        For each selected step, shows infected (top row) and recovered (bottom row) densities.
        Model parameters are displayed in a vertical panel on the right spanning both rows.
        """
        if step_indices is None:
            total_Rteps = len(self.I_files)
            step_indices = [0, 1 * total_Rteps // 4 - 1, 2 * total_Rteps // 4 - 1, total_Rteps - 1]
            step_indices = [max(0, min(idx, total_Rteps - 1)) for idx in step_indices]

        # Create figure with 2 rows and 5 columns (4 for snapshots + 1 for parameters)
        fig = plt.figure(figsize=(56, 16))  # taller to accommodate two rows
        gs = gridspec.GridSpec(2, 5, width_ratios=[1, 1, 1, 1, 1], wspace=0.3, hspace=0.1)
        for idx, step_idx in enumerate(step_indices[:4]):
            if step_idx >= len(self.I_files):
                continue

            # Load data
            I = self.load_matrix(self.I_files[step_idx])
            R = self.load_matrix(self.R_files[step_idx])

            # Infected subplot (row 0)
            ax_I = plt.subplot(gs[0, idx])
            im_I = ax_I.imshow(I, cmap='Reds', aspect='auto', origin='lower', vmin=0)
            ax_I.set_title(f'Step: {self.extract_step_number(self.I_files[step_idx])}', fontsize=36, fontweight='bold')
            #ax_I.set_xlabel('Antigenic coordinate (x)', fontsize=28)
            if step_idx == 0 :
                ax_I.set_ylabel('Antigenic coordinate (y)', fontsize=28)
            
            ax_I.tick_params(axis='both', labelsize=30)
            ax_I.set_facecolor('white')
            cbar_I = plt.colorbar(im_I, ax=ax_I)
            cbar_I.set_label('Infected density', fontsize=28)
            cbar_I.ax.tick_params(labelsize=24) 

            # Recovered subplot (row 1)
            ax_R = plt.subplot(gs[1, idx])
            im_R = ax_R.imshow(R, cmap='Greens', aspect='auto', origin='lower', vmin=0)
            #ax_R.set_title(f'Step: {self.extract_step_number(self.I_files[step_idx])}', fontsize=36, fontweight='bold')
            ax_R.set_xlabel('Antigenic coordinate (x)', fontsize=28)
            if step_idx == 0 : 
                ax_R.set_ylabel('Antigenic coordinate (y)', fontsize=28)

            ax_R.tick_params(axis='both', labelsize=30)
            ax_R.set_facecolor('white')
            cbar_R = plt.colorbar(im_R, ax=ax_R)
            cbar_R.set_label('Recovered density', fontsize=28)
            cbar_R.ax.tick_params(labelsize=24) 

        # Parameters panel on the right, spanning both rows
        ax_params = plt.subplot(gs[:, 4])   # all rows, last column
        ax_params.axis('off')
        params_text = self._format_parameters()
        ax_params.text(0.05, 0.5, params_text, transform=ax_params.transAxes, fontsize=30, fontfamily='monospace', verticalalignment='center', horizontalalignment='left', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='gray'))

        plt.suptitle('Virus wave evolution', fontsize=36, fontweight='bold')
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Snapshots with parameters saved as {save_path}")

    def normalize_matrix(self, matrix):
        """
        Приводит все элементы матрицы к диапазону [0, 1] на основе глобальных min и max.

        Параметры:
        matrix : array_like
            Входная матрица (любой тип, приводимый к float).

        Возвращает:
        numpy.ndarray
            Нормализованная матрица. Если все элементы исходной матрицы равны,
            возвращается нулевая матрица той же формы (чтобы избежать деления на ноль).
        """
        # Преобразуем входные данные в массив float (double)
        matrix = np.asarray(matrix, dtype=float)

        # Проверка на пустую матрицу
        if matrix.size == 0:
            return matrix

        min_val = matrix.min()
        max_val = matrix.max()

        # Если все значения одинаковы, нормализация не определена — возвращаем нули
        if max_val == min_val:
            return np.zeros_like(matrix)

        # Классическая минимаксная нормализация
        normalized = (matrix - min_val) / (max_val - min_val)
        return normalized

    def _format_parameters(self):
        """Formats model parameters into a multi-line string for display (in English)."""
        display_names = {
            'L': 'Grid size L',
            'M': 'Total steps M',
            'R0': 'Reproductive number R₀',
            'Ub': 'Mutation probability Ub',
            'a': 'Immunity distance a',
            'Varx': 'Coordinate spread Varx',
            'asym': 'Asymmetry type',
            'dt': 'Time step dt',
            'seed': 'Random seed',
        }

        order = [ 'L', 'M', 'R0', 'Ub', 'a', 'Varx', 'asym', 'dt', 'seed' ]

        lines = []
        for key in order:
            if key in self.model_params:
                value = self.model_params[key]
                if isinstance(value, (int, float)):
                    if isinstance(value, float):
                        if abs(value) < 1e-4 or abs(value) > 1e4:
                            display_value = f"{value:.2e}"
                        else:
                            # Remove trailing zeros after decimal
                            display_value = f"{value:.6f}".rstrip('0').rstrip('.') if '.' in f"{value:.6f}" else f"{value:.6f}"
                    else:
                        display_value = str(value)
                else:
                    display_value = str(value)

                name = display_names.get(key, key)
                lines.append(f"{name}: {display_value}")

        return "\n".join(lines)
    
    def create_wave_animation(self, save_path="wave_evolution.gif", max_frames=20):
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
        ax3 = plt.subplot(gs[1, 0])
        im3 = ax3.imshow(R0, cmap='Greens', aspect='auto', origin='lower', animated=True)
        ax3.set_title('Susceptible', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Antigenic coordinate (x)')
        ax3.set_ylabel('Antigenic coordinate (y)')
        plt.colorbar(im3, ax=ax3, label='Susceptible density')
        
        # 4. Model parameters (bottom right)
        ax4 = plt.subplot(gs[1, 1])
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
            composite = I * 0.8 + R * 0.4
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
            anim.save(save_path, writer=PillowWriter(fps=15))
        else:
            anim.save(save_path, writer='ffmpeg', fps=15)
        
        print(f"Animation saved as {save_path}")
        
        return anim

# Main function to run visualization
def main():
    parser = argparse.ArgumentParser(description='Visualization of antigenic evolution waves')
    parser.add_argument('-d', '--data-dir', type=str, default='output', help='Directory with data files (default: output)')
    parser.add_argument('-a', '--animation', type=str, default='wave_evolution.gif', help='Filename for animation output (default: wave_evolution.gif)')
    parser.add_argument('-s', '--snapshot', type=str, default='wave_snapshots.png', help='Filename for snapshot output (default: wave_snapshots.png)')
    parser.add_argument('--max-frames', type=int, default=3000, help='Maximum number of frames in animation (default: 3000)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("VIRUS ANTIGENIC EVOLUTION WAVE VISUALIZATION")
    print("=" * 60)
    print(f"Data directory: {args.data_dir}")
    print(f"Animation file: {args.animation}")
    print(f"Snapshot file: {args.snapshot}")
    print(f"Max frames: {args.max_frames}")
    
    visualizer = VirusWaveVisualizer(args.data_dir)
    
    if len(visualizer.I_files) == 0:
        print("No data to visualize!")
        return

    # 1. Create wave snapshots
    print("\n1. Creating wave snapshots...")
    try:
        visualizer.create_wave_snapshot(save_path=args.snapshot)
        print("✓ Snapshots created successfully")
    except Exception as e:
        print(f"✗ Error creating snapshots: {e}")
        import traceback
        traceback.print_exc()

    # 2. Create wave animation
    print("\n2. Creating wave animation...")
    try:
        anim = visualizer.create_wave_animation(save_path=args.animation, max_frames=args.max_frames)
        print("✓ Animation created successfully")
    except Exception as e:
        print(f"✗ Error creating animation: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 60)
    print("VISUALIZATION COMPLETE!")
    print("=" * 60)

if __name__ == "__main__":
    main()