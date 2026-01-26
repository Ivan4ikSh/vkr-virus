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
        
        if not os.path.exists(output_dir):
            print(f"Ошибка: Папка {output_dir} не существует!")
            return
            
        self.find_data_files()
        
        # Создаем кастомную цветовую карту для отображения "комет"
        colors = ['#FFFFFF', '#87CEEB', '#32CD32', '#FFD700', '#FF4500', '#8B0000']
        self.wave_cmap = LinearSegmentedColormap.from_list('wave', colors, N=256)
        
    def find_data_files(self):
        """Находит все файлы с результатами в папке output"""
        self.I_files = sorted(glob.glob(os.path.join(self.output_dir, "state_I_step_*.csv")), 
                            key=self.extract_step_number)
        self.R_files = sorted(glob.glob(os.path.join(self.output_dir, "state_R_step_*.csv")), 
                            key=self.extract_step_number)
        
        print(f"Найдено {len(self.I_files)} шагов симуляции")
        
    def extract_step_number(self, filename):
        """Извлекает номер шага из имени файла"""
        match = re.search(r'step_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    def load_matrix(self, filename):
        """Загружает матрицу из CSV файла"""
        return np.loadtxt(filename, delimiter=',')
    
    def create_wave_composite(self, I, R, step_num):
        """
        Создает композитное изображение волны (кометы)
        I - матрица инфицированных (голова кометы)
        R - матрица выздоровевших (хвост кометы)
        """
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # 1. Композитное изображение волны
        composite = np.zeros_like(I)
        # Голова (инфицированные) - высокие значения
        composite += I * 0.8
        # Хвост (выздоровевшие) - средние значения
        composite += R * 0.4
        
        im1 = axes[0].imshow(composite, cmap=self.wave_cmap, 
                            aspect='auto', origin='lower', 
                            interpolation='gaussian')
        axes[0].set_title(f'Волна эволюции (комета)\nШаг: {step_num}')
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
        
        axes[1].set_title('Контуры волны\nКрасный: инфицированные\nЗеленый: выздоровевшие')
        axes[1].set_xlabel('Антигенная координата (x) →')
        axes[1].set_ylabel('Нейтральная координата (y)')
        axes[1].grid(True, alpha=0.3)
        
        # Добавляем легенду для контуров
        axes[1].plot([], [], color='red', linewidth=2, label='Инфицированные (голова)')
        axes[1].plot([], [], color='green', linewidth=2, linestyle='--', label='Выздоровевшие (хвост)')
        axes[1].legend()
        
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
            
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Загружаем первый кадр
        I0 = self.load_matrix(self.I_files[0])
        R0 = self.load_matrix(self.R_files[0])
        
        # 1. Композитное изображение (полноцветная волна)
        composite0 = I0 * 0.8 + R0 * 0.4
        im1 = axes[0, 0].imshow(composite0, cmap=self.wave_cmap, 
                               aspect='auto', origin='lower', 
                               interpolation='gaussian', animated=True)
        axes[0, 0].set_title('Волна эволюции (комета)')
        axes[0, 0].set_xlabel('Антигенная координата (x)')
        axes[0, 0].set_ylabel('Нейтральная координата (y)')
        axes[0, 0].grid(True, alpha=0.3, linestyle='--')
        
        # 2. Инфицированные (голова)
        im2 = axes[0, 1].imshow(I0, cmap='Reds', aspect='auto', 
                               origin='lower', animated=True)
        axes[0, 1].set_title('Инфицированные (голова кометы)')
        axes[0, 1].set_xlabel('Антигенная координата (x)')
        axes[0, 1].set_ylabel('Нейтральная координата (y)')
        
        # 3. Выздоровевшие (хвост)
        im3 = axes[1, 0].imshow(R0, cmap='Greens', aspect='auto', 
                               origin='lower', animated=True)
        axes[1, 0].set_title('Выздоровевшие (хвост кометы)')
        axes[1, 0].set_xlabel('Антигенная координата (x)')
        axes[1, 0].set_ylabel('Нейтральная координата (y)')
        
        # 4. Логарифмическая разность (показатель "остроты" волны)
        diff0 = np.log1p(I0) - np.log1p(R0 + 1e-10)
        im4 = axes[1, 1].imshow(diff0, cmap='coolwarm', aspect='auto', 
                               origin='lower', animated=True)
        axes[1, 1].set_title('Градиент волны (лог-разность)')
        axes[1, 1].set_xlabel('Антигенная координата (x)')
        axes[1, 1].set_ylabel('Нейтральная координата (y)')
        
        plt.colorbar(im1, ax=axes[0, 0])
        plt.colorbar(im2, ax=axes[0, 1])
        plt.colorbar(im3, ax=axes[1, 0])
        plt.colorbar(im4, ax=axes[1, 1])
        
        fig.suptitle('Двумерная модель антигенной эволюции вируса\nШаг: 0', fontsize=14)
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
            fig.suptitle(f'Двумерная модель антигенной эволюции вируса\nШаг: {step_num}', fontsize=14)
            
            return images
        
        anim = FuncAnimation(fig, update, frames=min(20, len(self.I_files)), 
                           interval=200, blit=True)
        
        print(f"Создание анимации...")
        anim.save(save_path, writer=PillowWriter(fps=5))
        print(f"Анимация сохранена как {save_path}")
        
        return anim
    
    def analyze_wave_movement(self):
        """
        Анализирует движение волны и строит графики траекторий
        """
        if len(self.I_files) < 2:
            print("Недостаточно данных для анализа движения!")
            return
            
        # Находим центры масс инфицированных для каждого шага
        centers_x = []
        centers_y = []
        
        for I_file in self.I_files:
            I = self.load_matrix(I_file)
            
            # Вычисляем центр масс инфицированных
            total = I.sum()
            if total > 0:
                y_coords, x_coords = np.indices(I.shape)
                center_x = np.sum(x_coords * I) / total
                center_y = np.sum(y_coords * I) / total
            else:
                center_x, center_y = 0, 0
                
            centers_x.append(center_x)
            centers_y.append(center_y)
        
        # Строим график траектории
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Траектория центра масс
        axes[0].plot(centers_x, centers_y, 'b-', linewidth=2, alpha=0.7)
        axes[0].scatter(centers_x, centers_y, c=range(len(centers_x)), 
                       cmap='viridis', s=50, edgecolors='black')
        axes[0].set_title('Траектория центра масс волны')
        axes[0].set_xlabel('Антигенная координата (x)')
        axes[0].set_ylabel('Нейтральная координата (y)')
        axes[0].grid(True, alpha=0.3)
        axes[0].axis('equal')
        
        # Скорость движения волны
        steps = np.arange(len(centers_x))
        dx = np.diff(centers_x)
        dy = np.diff(centers_y)
        speed = np.sqrt(dx**2 + dy**2)
        
        axes[1].plot(steps[1:], speed, 'r-', linewidth=2)
        axes[1].set_title('Скорость движения волны')
        axes[1].set_xlabel('Шаг времени')
        axes[1].set_ylabel('Скорость (клеток/шаг)')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('output/wave_trajectory.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        return centers_x, centers_y
    
    def create_wave_snapshot(self, step_indices=None):
        """
        Создает снимки волны на разных этапах эволюции
        """
        if step_indices is None:
            step_indices = [0, len(self.I_files)//4, len(self.I_files)//2, 
                          3*len(self.I_files)//4, len(self.I_files)-1]
        
        n_steps = len(step_indices)
        fig, axes = plt.subplots(1, n_steps, figsize=(4*n_steps, 4))
        
        if n_steps == 1:
            axes = [axes]
        
        for idx, step_idx in enumerate(step_indices):
            if step_idx >= len(self.I_files):
                continue
                
            I = self.load_matrix(self.I_files[step_idx])
            R = self.load_matrix(self.R_files[step_idx])
            
            # Композитное изображение волны
            composite = I * 0.7 + R * 0.3
            
            im = axes[idx].imshow(composite, cmap=self.wave_cmap, 
                                 aspect='auto', origin='lower')
            step_num = self.extract_step_number(self.I_files[step_idx])
            axes[idx].set_title(f'Шаг {step_num}')
            axes[idx].set_xlabel('x')
            if idx == 0:
                axes[idx].set_ylabel('y')
        
        plt.suptitle('Эволюция волны (кометы) во времени', fontsize=14)
        plt.tight_layout()
        plt.savefig('output/wave_snapshots.png', dpi=150, bbox_inches='tight')
        plt.show()

# Основная функция для запуска визуализации
def main():
    print("=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ВОЛН АНТИГЕННОЙ ЭВОЛЮЦИИ ВИРУСА (КОМЕТЫ)")
    print("=" * 60)
    
    visualizer = VirusWaveVisualizer("output")
    
    if len(visualizer.I_files) == 0:
        print("Нет данных для визуализации!")
        return
    
    # 1. Анализ движения волны
    print("\n1. Анализ движения волны...")
    visualizer.analyze_wave_movement()
    
    # 2. Создаем анимацию волны
    print("\n2. Создание анимации волны...")
    try:
        visualizer.create_wave_animation("wave_evolution.gif")
    except Exception as e:
        print(f"Ошибка при создании анимации: {e}")
    
    # 3. Создаем снимки волны
    print("\n3. Создание снимков волны...")
    visualizer.create_wave_snapshot()

    print("\n" + "=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print("=" * 60)

if __name__ == "__main__":
    if not os.path.exists("output"):
        os.makedirs("output")
    
    main()