import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os
import glob
import re
import argparse
from matplotlib.colors import LinearSegmentedColormap, Normalize

class VirusWaveVisualizer:
    def __init__(self, data_dir="output"):
        """
        Инициализация визуализатора для модели волн вирусной эволюции
        
        Args:
            data_dir (str): Директория с файлами данных
        """
        self.data_dir = data_dir
        
        if not os.path.exists(data_dir):
            print(f"Ошибка: Папка {data_dir} не существует!")
            return
            
        self.find_data_files()
        
        # Создаем кастомную цветовую карту для отображения "комет"
        colors = ['#FFFFFF', '#87CEEB', '#32CD32', '#FFD700', '#FF4500', '#8B0000']
        self.wave_cmap = LinearSegmentedColormap.from_list('wave', colors, N=256)
        
    def find_data_files(self):
        """Находит все файлы с результатами в папке output"""
        self.I_files = sorted(glob.glob(os.path.join(self.data_dir, "state_I_step_*.csv")), key=self.extract_step_number)
        self.R_files = sorted(glob.glob(os.path.join(self.data_dir, "state_R_step_*.csv")), key=self.extract_step_number)
        
        if len(self.I_files) > 0:
            print(f"Найдено {len(self.I_files)} шагов симуляции в директории '{self.data_dir}'")
        else:
            print(f"Файлы данных не найдены в директории '{self.data_dir}'")
        
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
        
        im1 = axes[0].imshow(composite, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian')
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
        
        Args:
            save_path (str): Путь для сохранения анимации
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
        
        #anim = FuncAnimation(fig, update, frames=min(20, len(self.I_files)), interval=200, blit=True)
        
        print(f"Создание анимации...")
        anim.save(save_path, writer=PillowWriter(fps=5))
        print(f"Анимация сохранена как {save_path}")
        
        return anim
    
    
    def create_wave_snapshot(self, step_indices=None, save_path="wave_snapshots.png"):
        """
        Создает снимки волны на разных этапах эволюции
        
        Args:
            step_indices (list): Индексы шагов для создания снимков
            save_path (str): Путь для сохранения снимков
        """
        if step_indices is None:
            step_indices = [0, len(self.I_files)//4, len(self.I_files)//2, 3*len(self.I_files)//4, len(self.I_files)-1]
        
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
            
            im = axes[idx].imshow(composite, cmap=self.wave_cmap, aspect='auto', origin='lower')
            step_num = self.extract_step_number(self.I_files[step_idx])
            axes[idx].set_title(f'Шаг {step_num}')
            axes[idx].set_xlabel('x')
            if idx == 0:
                axes[idx].set_ylabel('y')
        
        plt.suptitle(f'Эволюция волны (кометы) во времени\nДиректория: {self.data_dir}', fontsize=14)
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        #plt.show()

# Основная функция для запуска визуализации
def main():
    parser = argparse.ArgumentParser(description='Визуализация волн антигенной эволюции вируса')
    parser.add_argument('-d', '--data-dir', type=str, default='output', help='Директория с файлами данных (по умолчанию: output)')
    parser.add_argument('-a', '--animation', type=str, default='wave_evolution.gif', help='Имя файла для сохранения анимации (по умолчанию: wave_evolution.gif)')
    parser.add_argument('-s', '--snapshot', type=str, default='wave_snapshots.png', help='Имя файла для сохранения снимков (по умолчанию: wave_snapshots.png)')
    parser.add_argument('--max-frames', type=int, default=20, help='Максимальное количество кадров в анимации (по умолчанию: 20)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ВОЛН АНТИГЕННОЙ ЭВОЛЮЦИИ ВИРУСА")
    print("=" * 60)
    print(f"Директория с данными: {args.data_dir}")
    print(f"Файл анимации: {args.animation}")
    print(f"Файл снимков: {args.snapshot}")
    print(f"Максимальное количество кадров: {args.max_frames}")
    
    visualizer = VirusWaveVisualizer(args.data_dir)
    
    if len(visualizer.I_files) == 0:
        print("Нет данных для визуализации!")
        return

    # 1. Создаем анимацию волны
    print("\n1 Создание анимации волны...")
    try:
        # Временно модифицируем метод create_wave_animation для использования max_frames
        # Создаем временную функцию, или модифицируем класс
        # В данном случае, проще всего создать анимацию с ограничением кадров
        
        # Загружаем данные для анимации
        if len(visualizer.I_files) > 0:
            # Используем тот же код, что и в create_wave_animation, но с параметром max_frames
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            
            I0 = visualizer.load_matrix(visualizer.I_files[0])
            R0 = visualizer.load_matrix(visualizer.R_files[0])
            
            composite0 = I0 * 0.8 + R0 * 0.4
            im1 = axes[0, 0].imshow(composite0, cmap=visualizer.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian', animated=True)
            axes[0, 0].set_title('Волна эволюции (комета)')
            axes[0, 0].set_xlabel('Антигенная координата (x)')
            axes[0, 0].set_ylabel('Нейтральная координата (y)')
            axes[0, 0].grid(True, alpha=0.3, linestyle='--')
            
            im2 = axes[0, 1].imshow(I0, cmap='Reds', aspect='auto', origin='lower', animated=True)
            axes[0, 1].set_title('Инфицированные (голова кометы)')
            axes[0, 1].set_xlabel('Антигенная координата (x)')
            axes[0, 1].set_ylabel('Нейтральная координата (y)')
            
            im3 = axes[1, 0].imshow(R0, cmap='Greens', aspect='auto', origin='lower', animated=True)
            axes[1, 0].set_title('Выздоровевшие (хвост кометы)')
            axes[1, 0].set_xlabel('Антигенная координата (x)')
            axes[1, 0].set_ylabel('Нейтральная координата (y)')
            
            diff0 = np.log1p(I0) - np.log1p(R0 + 1e-10)
            im4 = axes[1, 1].imshow(diff0, cmap='coolwarm', aspect='auto', origin='lower', animated=True)
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
                I = visualizer.load_matrix(visualizer.I_files[frame])
                R = visualizer.load_matrix(visualizer.R_files[frame])
                
                composite = I * 0.8 + R * 0.4
                images[0].set_data(composite)
                images[1].set_data(I)
                images[2].set_data(R)
                
                diff = np.log1p(I) - np.log1p(R + 1e-10)
                images[3].set_data(diff)
                
                for im in images:
                    im.autoscale()
                
                step_num = visualizer.extract_step_number(visualizer.I_files[frame])
                fig.suptitle(f'Двумерная модель антигенной эволюции вируса\nШаг: {step_num}', fontsize=14)
                
                return images
            
            max_frames = min(args.max_frames, len(visualizer.I_files))
            anim = FuncAnimation(fig, update, frames=max_frames, interval=200, blit=True)
            
            print(f"Создание анимации с {max_frames} кадрами...")
            anim.save(args.animation, writer=PillowWriter(fps=5))
            print(f"Анимация сохранена как {args.animation}")
    except Exception as e:
        print(f"Ошибка при создании анимации: {e}")
    
    # 2. Создаем снимки волны
    print("\n2 Создание снимков волны...")
    try:
        visualizer.create_wave_snapshot(save_path=args.snapshot)
    except Exception as e:
        print(f"Ошибка при создании снимков: {e}")

    print("\n" + "=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print("=" * 60)

if __name__ == "__main__":
    main()