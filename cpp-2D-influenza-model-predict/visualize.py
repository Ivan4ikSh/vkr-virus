import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob
import re
from matplotlib.colors import LinearSegmentedColormap, Normalize

class VirusWaveVisualizer:
    def __init__(self, output_dir="output"):
        """
        Инициализация визуализатора для модели волн вирусной эволюции
        """
        self.output_dir = output_dir
        self.simulation_params = {}
        
        if not os.path.exists(output_dir):
            print(f"Ошибка: Папка {output_dir} не существует!")
            return
            
        self.find_data_files()
        self.load_simulation_parameters()
        
        # Создаем кастомную цветовую карту для отображения "комет"
        colors = ['#FFFFFF', '#87CEEB', '#32CD32', '#FFD700', '#FF4500', '#8B0000']
        self.wave_cmap = LinearSegmentedColormap.from_list('wave', colors, N=256)
        
    def find_data_files(self):
        """Находит все файлы с результатами в папке output"""
        self.I_files = sorted(glob.glob(os.path.join(self.output_dir, "state_I_step_*.csv")), key=self.extract_step_number)
        self.R_files = sorted(glob.glob(os.path.join(self.output_dir, "state_R_step_*.csv")), key=self.extract_step_number)
        
        print(f"Найдено {len(self.I_files)} шагов симуляции")
        
    def extract_step_number(self, filename):
        """Извлекает номер шага из имени файла"""
        match = re.search(r'step_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    def load_matrix(self, filename):
        """Загружает матрицу из CSV файла"""
        return np.loadtxt(filename, delimiter=',')
    
    def load_simulation_parameters(self):
        """Загружает параметры симуляции из файла или устанавливает значения по умолчанию"""
        params_file = os.path.join(self.output_dir, "simulation_params.txt")
        
        if os.path.exists(params_file):
            print("Загрузка параметров симуляции...")
            with open(params_file, 'r') as f:
                for line in f:
                    if '=' in line:
                        key, value = line.strip().split('=', 1)
                        self.simulation_params[key.strip()] = value.strip()
        else:
            # Если файла нет, используем значения по умолчанию
            print("Файл параметров не найден, используются значения по умолчанию")
            self.simulation_params = {
                'R0': '2.6',
                'D': '0.0001',
                'a': '7',
                'N': '1e10',
                'Varx': '0.01',
                'L': '50',
                'Tmax': '600',
                'stept': '1',
                'asymmetry': 'ANI_ASY'
            }
            
            # Сохраняем параметры по умолчанию в файл
            self.save_default_parameters()
    
    def save_default_parameters(self):
        """Сохраняет параметры по умолчанию в файл"""
        params_file = os.path.join(self.output_dir, "simulation_params.txt")
        with open(params_file, 'w') as f:
            for key, value in self.simulation_params.items():
                f.write(f"{key}={value}\n")
        print(f"Параметры сохранены в {params_file}")
    
    def get_params_string(self, include_all=False):
        """Возвращает строку с параметрами симуляции для заголовка"""
        if include_all:
            # Полный набор параметров
            params = [
                f"R0={self.simulation_params.get('R0', '?')}",
                f"D={self.simulation_params.get('D', '?')}",
                f"a={self.simulation_params.get('a', '?')}",
                f"N={self.simulation_params.get('N', '?')}",
                f"Varx={self.simulation_params.get('Varx', '?')}",
                f"L={self.simulation_params.get('L', '?')}",
                f"Tmax={self.simulation_params.get('Tmax', '?')}"
            ]
        else:
            # Основные параметры
            params = [
                f"R0={self.simulation_params.get('R0', '?')}",
                f"D={self.simulation_params.get('D', '?')}",
                f"a={self.simulation_params.get('a', '?')}",
                f"N={self.simulation_params.get('N', '?')}"
            ]
        
        return ", ".join(params)
    
    def create_wave_composite(self, I, R, step_num):
        """
        Создает композитное изображение волны (кометы)
        I - матрица инфицированных (голова кометы)
        R - матрица выздоровевших (хвост кометы)
        """
        params_str = self.get_params_string()
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        
        # Общий заголовок с параметрами
        main_title = f"Двумерная модель антигенной эволюции вируса\nШаг: {step_num}, Параметры: {params_str}"
        fig.suptitle(main_title, fontsize=14, y=1.02)
        
        # 1. Композитное изображение волны
        composite = np.zeros_like(I)
        # Голова (инфицированные) - высокие значения
        composite += I * 0.8
        # Хвост (выздоровевшие) - средние значения
        composite += R * 0.4
        
        im1 = axes[0].imshow(composite, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian')
        axes[0].set_title('Волна эволюции (комета)')
        axes[0].set_xlabel('Антигенная координата (x) →')
        axes[0].set_ylabel('Нейтральная координата (y)')
        axes[0].grid(True, alpha=0.3, linestyle='--')
        plt.colorbar(im1, ax=axes[0], label='Интенсивность волны')
        
        # 2. Контурное представление волны
        X, Y = np.meshgrid(np.arange(I.shape[1]), np.arange(I.shape[0]))
        
        # Контуры инфицированных (голова)
        contours_I = axes[1].contour(X, Y, I, levels=5, colors='red', linewidths=2, alpha=0.7)
        # Контуры выздоровевших (хвост)
        contours_R = axes[1].contour(X, Y, R, levels=5, colors='green', linewidths=2, alpha=0.7, linestyles='--')
        
        axes[1].set_title('Контуры волны')
        axes[1].set_xlabel('Антигенная координата (x) →')
        axes[1].set_ylabel('Нейтральная координата (y)')
        axes[1].grid(True, alpha=0.3)
        
        # Добавляем легенду для контуров
        axes[1].plot([], [], color='red', linewidth=2, label='Инфицированные (голова)')
        axes[1].plot([], [], color='green', linewidth=2, linestyle='--', label='Выздоровевшие (хвост)')
        axes[1].legend(loc='upper right', fontsize=9)
        
        # 3. 3D представление волны
        ax3d = fig.add_subplot(133, projection='3d')
        
        # Нормализуем данные для лучшей визуализации
        I_norm = I / (I.max() + 1e-10)
        R_norm = R / (R.max() + 1e-10)
        
        surf1 = ax3d.plot_surface(X, Y, I_norm, cmap='Reds', alpha=0.7, label='Инфицированные')
        surf2 = ax3d.plot_surface(X, Y, R_norm, cmap='Greens', alpha=0.5, label='Выздоровевшие')
        
        ax3d.set_title('3D профиль волны')
        ax3d.set_xlabel('Антигенная координата (x)')
        ax3d.set_ylabel('Нейтральная координата (y)')
        ax3d.set_zlabel('Плотность')
        
        plt.tight_layout()
        return fig, composite
    
    def create_wave_animation(self, save_path="wave_evolution.gif"):
        """
        Создает анимацию эволюции волны (кометы)
        """
        if len(self.I_files) == 0:
            print("Нет данных для анимации!")
            return
            
        params_str = self.get_params_string()
        
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        
        # Общий заголовок с параметрами
        fig.suptitle(f'Двумерная модель антигенной эволюции вируса\nПараметры: {params_str}', fontsize=14, y=1.02)
        
        # Загружаем первый кадр
        I0 = self.load_matrix(self.I_files[0])
        R0 = self.load_matrix(self.R_files[0])
        
        # 1. Композитное изображение (полноцветная волна)
        composite0 = I0 * 0.8 + R0 * 0.4
        im1 = axes[0, 0].imshow(composite0, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian', animated=True)
        axes[0, 0].set_title('Волна эволюции (комета)')
        axes[0, 0].set_xlabel('Антигенная координата (x)')
        axes[0, 0].set_ylabel('Нейтральная координата (y)')
        axes[0, 0].grid(True, alpha=0.3, linestyle='--')
        
        # 2. Инфицированные (голова)
        im2 = axes[0, 1].imshow(I0, cmap='Reds', aspect='auto', origin='lower', animated=True)
        axes[0, 1].set_title('Инфицированные (голова кометы)')
        axes[0, 1].set_xlabel('Антигенная координата (x)')
        axes[0, 1].set_ylabel('Нейтральная координата (y)')
        
        # 3. Выздоровевшие (хвост)
        im3 = axes[1, 0].imshow(R0, cmap='Greens', aspect='auto', origin='lower', animated=True)
        axes[1, 0].set_title('Выздоровевшие (хвост кометы)')
        axes[1, 0].set_xlabel('Антигенная координата (x)')
        axes[1, 0].set_ylabel('Нейтральная координата (y)')
        
        # 4. Логарифмическая разность (показатель "остроты" волны)
        diff0 = np.log1p(I0) - np.log1p(R0 + 1e-10)
        im4 = axes[1, 1].imshow(diff0, cmap='coolwarm', aspect='auto', origin='lower', animated=True)
        axes[1, 1].set_title('Градиент волны (лог-разность)')
        axes[1, 1].set_xlabel('Антигенная координата (x)')
        axes[1, 1].set_ylabel('Нейтральная координата (y)')
        
        plt.colorbar(im1, ax=axes[0, 0], label='Интенсивность')
        plt.colorbar(im2, ax=axes[0, 1], label='Плотность')
        plt.colorbar(im3, ax=axes[1, 0], label='Плотность')
        plt.colorbar(im4, ax=axes[1, 1], label='Градиент')
        
        plt.tight_layout()
        
        images = [im1, im2, im3, im4]
        
        def update(frame):
            I = self.load_matrix(self.I_files[frame])
            R = self.load_matrix(self.R_files[frame])
            
            # Обновляем все изображения
            composite = I * 0.8 + R * 0.4
            images[0].set_data(composite)
            images[1].set_data(I)
            images[2].set_data(R)
            
            # Логарифмическая разность
            diff = np.log1p(I) - np.log1p(R + 1e-10)
            images[3].set_data(diff)
            
            # Автомасштабирование
            for im in images:
                im.autoscale()
            
            step_num = self.extract_step_number(self.I_files[frame])
            # Обновляем заголовок для текущего шага
            fig.suptitle(f'Двумерная модель антигенной эволюции вируса\nШаг: {step_num}, Параметры: {params_str}', fontsize=14, y=1.02)
            
            return images
        
        num_frames = min(30, len(self.I_files))  # Ограничиваем количество кадров для скорости
        anim = FuncAnimation(fig, update, frames=num_frames, interval=300, blit=True)
        
        print(f"Создание анимации ({num_frames} кадров)...")
        
        # Сохраняем анимацию
        try:
            anim.save(save_path, writer=PillowWriter(fps=5))
            print(f"Анимация сохранена как {save_path}")
        except Exception as e:
            print(f"Ошибка при сохранении анимации: {e}")
            # Альтернативно показываем анимацию в окне
            plt.show()
        
        return anim
    
    def create_wave_snapshot(self, step_indices=None):
        """
        Создает снимки волны на разных этапах эволюции
        """
        if step_indices is None:
            step_indices = [0, len(self.I_files)//4, len(self.I_files)//2, 3*len(self.I_files)//4, len(self.I_files)-1]
        
        n_steps = len(step_indices)
        fig, axes = plt.subplots(1, n_steps, figsize=(4*n_steps, 4))
        
        if n_steps == 1:
            axes = [axes]
        
        params_str = self.get_params_string(include_all=True)
        
        for idx, step_idx in enumerate(step_indices):
            if step_idx >= len(self.I_files):
                continue
                
            I = self.load_matrix(self.I_files[step_idx])
            R = self.load_matrix(self.R_files[step_idx])
            
            # Композитное изображение волны
            composite = I * 0.7 + R * 0.3
            
            im = axes[idx].imshow(composite, cmap=self.wave_cmap, aspect='auto', origin='lower')
            step_num = self.extract_step_number(self.I_files[step_idx])
            axes[idx].set_title(f'Шаг {step_num}')
            axes[idx].set_xlabel('Антигенная координата (x)')
            if idx == 0:
                axes[idx].set_ylabel('Нейтральная координата (y)')
        
        # Общий заголовок с параметрами
        plt.suptitle(f'Эволюция волны (кометы) во времени\nПараметры: {params_str}', fontsize=14, y=1.05)
        plt.tight_layout()
        
        output_path = 'output/wave_snapshots.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Снимки сохранены как {output_path}")
        plt.show()

# Основная функция для запуска визуализации
def main():
    visualizer = VirusWaveVisualizer("output")
    
    if len(visualizer.I_files) == 0:
        print("Нет данных для визуализации!")
        return

    # Выводим параметры симуляции
    print("\nПараметры симуляции:")
    for key, value in visualizer.simulation_params.items():
        print(f"  {key}: {value}")
    
    # 1. Создаем анимацию волны
    print("\n1. Создание анимации волны...")
    try:
        visualizer.create_wave_animation("wave_evolution.gif")
    except Exception as e:
        print(f"Ошибка при создании анимации: {e}")
    
    # 2. Создаем снимки волны
    print("\n2. Создание снимков волны...")
    visualizer.create_wave_snapshot()
    
    # 3. Создаем отдельные композитные изображения для первых 3 шагов
    print("\n3. Создание отдельных композитных изображений...")
    for i in range(min(3, len(visualizer.I_files))):
        I = visualizer.load_matrix(visualizer.I_files[i])
        R = visualizer.load_matrix(visualizer.R_files[i])
        step_num = visualizer.extract_step_number(visualizer.I_files[i])
        
        fig, composite = visualizer.create_wave_composite(I, R, step_num)
        
        output_path = f'output/wave_composite_step_{step_num}.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Сохранено: {output_path}")
        plt.close(fig)
        
        # Сохраняем композитное изображение
        plt.imsave(f'output/composite_image_step_{step_num}.png', composite, cmap=visualizer.wave_cmap)
        print(f"  Сохранено изображение: output/composite_image_step_{step_num}.png")

    print("\n" + "=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print("=" * 60)

if __name__ == "__main__":
    if not os.path.exists("output"):
        os.makedirs("output")
    
    main()