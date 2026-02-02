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
        
        # Загружаем параметры модели из файла
        self.model_params = self.load_model_parameters()
        
    def load_model_parameters(self):
        """Загружает параметры модели из файла model_parameters.txt"""
        param_file = os.path.join(self.data_dir, "model_parameters.txt")
        
        if not os.path.exists(param_file):
            print(f"Файл параметров не найден: {param_file}")
            # Возвращаем параметры по умолчанию
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
                    
                    # Парсим строки вида "Параметр: значение"
                    if ': ' in line:
                        key, value = line.split(': ', 1)
                        
                        # Преобразуем числовые значения
                        if key in ['Grid size L', 'Reproductive number R0', 
                                  'Mutation probability D', 'Immunity distance a',
                                  'Coordinate spread Varx', 'Initial infected fraction',
                                  'Show interval tshow', 'Start time T0', 'Total steps M']:
                            try:
                                # Обрабатываем научную нотацию
                                if 'e' in value.lower() or 'E' in value.lower():
                                    params[key] = float(value)
                                else:
                                    # Проверяем, является ли значение целым числом
                                    if '.' not in value and value.replace('.', '').isdigit():
                                        params[key] = int(value)
                                    else:
                                        params[key] = float(value)
                            except ValueError:
                                params[key] = value
                        else:
                            params[key] = value
            
            # Добавляем дополнительную информацию
            if len(self.I_files) > 0:
                params['Total simulation files'] = len(self.I_files)
                if len(self.I_files) > 0:
                    params['First step'] = self.extract_step_number(self.I_files[0])
                    params['Last step'] = self.extract_step_number(self.I_files[-1])
                    
                    # Загружаем один файл для определения размера
                    sample_data = self.load_matrix(self.I_files[0])
                    params['Data dimensions'] = f"{sample_data.shape[0]}x{sample_data.shape[1]}"
            
            print(f"✓ Загружено {len(params)} параметров модели")
            return params
            
        except Exception as e:
            print(f"Ошибка при загрузке параметров: {e}")
            return self.get_default_parameters()
    
    def get_default_parameters(self):
        """Возвращает параметры по умолчанию если файл не найден"""
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
    
    def create_wave_snapshot(self, step_indices=None, save_path="wave_snapshots.png"):
        """
        Создает снимки волны на разных этапах эволюции
        6 позиций: 5 графиков + таблица параметров
        """
        if step_indices is None:
            # Выбираем 5 равномерно распределенных шагов
            total_steps = len(self.I_files)
            step_indices = [0, max(1, total_steps // 4), total_steps // 2, 3 * total_steps // 4, max(0, total_steps - 1)]
        
        # Создаем фигуру с 6 позициями
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()  # Преобразуем в плоский массив
        
        # Обрабатываем 5 графиков
        for idx, step_idx in enumerate(step_indices[:5]):  # Берем только первые 5
            if step_idx >= len(self.I_files):
                continue
                
            I = self.load_matrix(self.I_files[step_idx])
            R = self.load_matrix(self.R_files[step_idx])
            
            # Композитное изображение волны
            composite = I * 0.8 + R * 0.4
            # Основной график
            ax = axes[idx]
            im = ax.imshow(composite, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian')
            ax.set_title(f'Шаг: {self.extract_step_number(self.I_files[step_idx])}', fontsize=14, fontweight='bold')
            ax.set_xlabel('Антигенная координата (x)', fontsize=10)
            ax.set_ylabel('Антигенная координата (y)', fontsize=10)
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Добавляем colorbar к каждому графику
            plt.colorbar(im, ax=ax, label='Интенсивность волны')
        
        # 6-я позиция - таблица с параметрами
        ax_params = axes[5]
        ax_params.axis('off')
        
        # Подготовка текста параметров
        params_text = "ПАРАМЕТРЫ МОДЕЛИ:\n\n"
        # Форматируем параметры
        param_display_order = ["Grid size L", "Reproductive number R0",  "Mutation probability D", "Immunity distance a", "Coordinate spread Varx", "Total population N", "Initial infected fraction", "Show interval tshow", "Start time T0",  "Total steps M", "Asymmetry type", "Total simulation files", "First step", "Last step", "Data dimensions"]
        
        for key in param_display_order:
            if key in self.model_params:
                value = self.model_params[key]
                # Форматируем значения
                if isinstance(value, float):
                    if abs(value) < 1e-4 or abs(value) > 1e4:
                        display_value = f"{value:.2e}"
                    else:
                        display_value = f"{value:.6f}"
                else:
                    display_value = str(value)
                
                # Преобразуем названия параметров на русский
                russian_names = {
                    "Grid size L": "Размер сетки",
                    "Reproductive number R0": "Число R₀",
                    "Mutation probability D": "Вероятность мутации D",
                    "Immunity distance a": "Иммунное расстояние a",
                    "Coordinate spread Varx": "Разброс координат Varx",
                    "Total population N": "Популяция N",
                    "Initial infected fraction": "Начальные инфицированные",
                    "Show interval tshow": "Интервал сохранения",
                    "Start time T0": "Начало записи T₀",
                    "Total steps M": "Всего шагов",
                    "Asymmetry type": "Тип асимметрии",
                    "Total simulation files": "Файлов сохранено",
                    "First step": "Первый шаг",
                    "Last step": "Последний шаг",
                    "Data dimensions": "Размер данных"
                }
                
                russian_key = russian_names.get(key, key)
                params_text += f"{russian_key}: {display_value}\n"
        
        # Добавляем текст с увеличенным шрифтом
        ax_params.text(0.1, 0.95, params_text, fontsize=12, fontfamily='monospace', verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9, edgecolor='navy', linewidth=2))
        # Добавляем общий заголовок
        plt.suptitle('Эволюция вирусной волны', fontsize=18, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Снимки сохранены как {save_path}")
        #plt.show()
    
    def create_wave_animation(self, save_path="wave_evolution.gif", max_frames=20):
        """
        Создает анимацию эволюции волны с 3 графиками и параметрами
        """
        if len(self.I_files) == 0:
            print("Нет данных для анимации!")
            return
            
        # Ограничиваем количество кадров
        total_frames = min(max_frames, len(self.I_files))
        
        fig = plt.figure(figsize=(18, 10))
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
        
        # Загружаем первый кадр для инициализации
        I0 = self.load_matrix(self.I_files[0])
        R0 = self.load_matrix(self.R_files[0])
        
        # 1. Композитное изображение (верхний левый)
        ax1 = plt.subplot(gs[0, 0])
        composite0 = I0 * 0.8 + R0 * 0.4
        im1 = ax1.imshow(composite0, cmap=self.wave_cmap, aspect='auto', origin='lower', interpolation='gaussian', animated=True)
        ax1.set_title('Волна эволюции (комета)', fontsize=12, fontweight='bold')
        ax1.set_xlabel('Антигенная координата (x)')
        ax1.set_ylabel('Нейтральная координата (y)')
        ax1.grid(True, alpha=0.3, linestyle='--')
        plt.colorbar(im1, ax=ax1, label='Интенсивность')
        
        # 2. Инфицированные (верхний правый)
        ax2 = plt.subplot(gs[0, 1])
        im2 = ax2.imshow(I0, cmap='Reds', aspect='auto', origin='lower', animated=True)
        ax2.set_title('Инфицированные', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Антигенная координата (x)')
        ax2.set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im2, ax=ax2, label='Плотность инфицированных')
        
        # 3. Выздоровевшие (нижний левый)
        ax3 = plt.subplot(gs[1, 0])
        im3 = ax3.imshow(R0, cmap='Greens', aspect='auto', origin='lower', animated=True)
        ax3.set_title('Выздоровевшие', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Антигенная координата (x)')
        ax3.set_ylabel('Нейтральная координата (y)')
        plt.colorbar(im3, ax=ax3, label='Плотность выздоровевших')
        
        # 4. Параметры модели (нижний правый)
        ax4 = plt.subplot(gs[1, 1])
        ax4.axis('off')
        
        # Подготовка текста параметров
        params_text = "ПАРАМЕТРЫ МОДЕЛИ:\n\n"
        
        # Отображаем только основные параметры
        main_params = [
            "Grid size L",
            "Reproductive number R0", 
            "Mutation probability D",
            "Immunity distance a",
            "Total population N",
            "Asymmetry type"
        ]
        
        russian_names = {
            "Grid size L": "Размер сетки",
            "Reproductive number R0": "Число R₀",
            "Mutation probability D": "Вероятность мутации D",
            "Immunity distance a": "Иммунное расстояние a",
            "Total population N": "Популяция N",
            "Asymmetry type": "Тип асимметрии"
        }
        
        for key in main_params:
            if key in self.model_params:
                value = self.model_params[key]
                if isinstance(value, float):
                    if abs(value) < 1e-4 or abs(value) > 1e4:
                        display_value = f"{value:.2e}"
                    else:
                        display_value = f"{value:.6f}"
                else:
                    display_value = str(value)
                
                russian_key = russian_names.get(key, key)
                params_text += f"{russian_key}: {display_value}\n"
        
        # Добавляем информацию о файлах и размере
        params_text += f"\nВсего шагов: {len(self.I_files)}\n"
        params_text += f"Размер данных: {I0.shape[0]}x{I0.shape[1]}"
        self.params_text_obj = ax4.text(0.05, 0.95, params_text, fontsize=15, fontfamily='monospace',  verticalalignment='top', transform=ax4.transAxes, bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        step_num = self.extract_step_number(self.I_files[0])
        images = [im1, im2, im3]
        
        def update(frame):
            # Загружаем данные для текущего кадра
            I = self.load_matrix(self.I_files[frame])
            R = self.load_matrix(self.R_files[frame])
            
            # Обновляем изображения
            composite = I * 0.8 + R * 0.4
            images[0].set_data(composite)
            images[1].set_data(I)
            images[2].set_data(R)
            
            # Автомасштабирование
            for im in images:
                im.autoscale()
            
            # Обновляем информацию о шаге
            step_num = self.extract_step_number(self.I_files[frame])            
            # Обновляем заголовок
            fig.suptitle(f'Анимация эволюции вирусной волны - Шаг: {step_num}', fontsize=18, fontweight='bold')
            
            return images
        
        # Создаем анимацию
        anim = FuncAnimation(fig, update, frames=total_frames, interval=300, blit=True, repeat=True)
        
        print(f"Создание анимации с {total_frames} кадрами...")
        
        # Сохраняем анимацию
        if save_path.endswith('.gif'):
            anim.save(save_path, writer=PillowWriter(fps=15))
        else:
            anim.save(save_path, writer='ffmpeg', fps=15)
        
        print(f"Анимация сохранена как {save_path}")
        
        return anim

# Основная функция для запуска визуализации
def main():
    parser = argparse.ArgumentParser(description='Визуализация волн антигенной эволюции вируса')
    parser.add_argument('-d', '--data-dir', type=str, default='output', help='Директория с файлами данных (по умолчанию: output)')
    parser.add_argument('-a', '--animation', type=str, default='wave_evolution.gif', help='Имя файла для сохранения анимации (по умолчанию: wave_evolution.gif)')
    parser.add_argument('-s', '--snapshot', type=str, default='wave_snapshots.png', help='Имя файла для сохранения снимков (по умолчанию: wave_snapshots.png)')
    parser.add_argument('--max-frames', type=int, default=3000, help='Максимальное количество кадров в анимации (по умолчанию: 300)')
    
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

    # 1. Создаем снимки волны
    print("\n1. Создание снимков волны...")
    try:
        visualizer.create_wave_snapshot(save_path=args.snapshot)
        print("✓ Снимки успешно созданы")
    except Exception as e:
        print(f"✗ Ошибка при создании снимков: {e}")
        import traceback
        traceback.print_exc()

    # 2. Создаем анимацию волны
    print("\n2. Создание анимации волны...")
    try:
        anim = visualizer.create_wave_animation(save_path=args.animation, max_frames=args.max_frames)
        print("✓ Анимация успешно создана")
    except Exception as e:
        print(f"✗ Ошибка при создании анимации: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 60)
    print("ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print("=" * 60)

if __name__ == "__main__":
    main()