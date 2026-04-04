import os
import sys
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import re

# Increase recursion limit for deep trees
sys.setrecursionlimit(20000)

def get_params_string(data_dir):
    """
    Reads parameters from observables.txt and formats them into a 2-line string for titles.
    """
    obs_file = os.path.join(data_dir, 'observables.txt')
    if not os.path.exists(obs_file):
        return ""
    
    params = []
    with open(obs_file, 'r') as f:
        for line in f:
            if '=' in line:
                params.append(line.strip())
    
    if not params:
        return ""
    
    # Разделяем параметры пополам для двух строк
    mid = len(params) // 2 + 1 if len(params) > 5 else len(params)
    line1 = ", ".join(params[:mid])
    line2 = ", ".join(params[mid:])
    
    if line2:
        return f"{line1}\n{line2}"
    return line1

def plot_timeseries(data_dir, output_dir, params_str):
    ts_file = os.path.join(data_dir, 'timeseries.csv')
    if not os.path.exists(ts_file):
        print(f"[Warning] File not found: {ts_file}")
        return

    df = pd.read_csv(ts_file)
    
    plt.figure(figsize=(10, 6))
    plt.plot(df['step'], df['total_i'], color='red', lw=2)
    
    # Добавляем параметры в заголовок
    title = 'Epidemic Dynamics: Total Fraction of Infected (finf)'
    if params_str: title += f'\n{params_str}'
    plt.title(title, fontsize=11)
    
    plt.xlabel('Simulation Step')
    plt.ylabel('Fraction of Infected (total_i)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    out_file = os.path.join(output_dir, 'timeseries.png')
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Timeseries plot saved: {out_file}")

def plot_wave(data_dir, output_dir, params_str):
    wave_file = os.path.join(data_dir, 'wave.csv')
    if not os.path.exists(wave_file):
        print(f"[Warning] File not found: {wave_file}")
        return

    df = pd.read_csv(wave_file)
    
    plt.figure(figsize=(14, 6))
    
    depths = df['depth'].unique()
    cmap = plt.get_cmap('tab20')
    
    for idx, d in enumerate(sorted(depths)):
        d_data = df[df['depth'] == d]
        plt.plot(d_data['step'], d_data['i'], color=cmap(idx % 20), lw=1.5)
        
    plt.yscale('log')
    plt.ylim(1e-8, 1e-1) 
    
    # Добавляем параметры в заголовок
    title = 'Evolution of Viral Lineages (Traveling Wave)'
    if params_str: title += f'\n{params_str}'
    plt.title(title, fontsize=12)
    
    plt.xlabel('Simulation Step')
    plt.ylabel('Fraction of Infected (log)')
    
    plt.grid(True, which="both", linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    out_file = os.path.join(output_dir, 'traveling_wave.png')
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Traveling wave plot saved: {out_file}")

def calculate_concentrated_layout(df):
    children = {}
    is_infected_dict = {}
    
    for _, row in df.iterrows():
        nid = int(row['id'])
        is_infected_dict[nid] = int(row['is_infected']) == 1
        children[nid] = []
        
    root_id = None
    for _, row in df.iterrows():
        pid = int(row['parent_id'])
        cid = int(row['id'])
        if pid == -1:
            root_id = cid
        elif pid in children:
            children[pid].append(cid)
            
    if root_id is None:
        root_id = 0 

    node_weights = {}

    def calculate_weights(node):
        weight = 1 if is_infected_dict[node] else 0
        for child in children[node]:
            weight += calculate_weights(child)
        node_weights[node] = weight
        return weight

    calculate_weights(root_id)

    node_y = {}
    current_y = [0.0]

    def dfs_layout(node):
        if not children.get(node):
            node_y[node] = current_y[0]
            current_y[0] += 1.0
        else:
            sorted_children = sorted(children[node], key=lambda c: node_weights[c])
            y_sum = 0.0
            for child in sorted_children:
                dfs_layout(child)
                y_sum += node_y[child]
            node_y[node] = y_sum / len(sorted_children)

    dfs_layout(root_id)
    return node_y

def find_main_trunk(df):
    if df.empty:
        return set(), set()
    
    # 1. Находим самый глубокий узел
    farthest_node_row = df.loc[df['depth'].idxmax()]
    farthest_node_id = int(farthest_node_row['id'])
    
    parent_of = {}
    is_infected_dict = {}
    children_of = {}
    
    # Считываем связи и статус инфекции
    for _, row in df.iterrows():
        nid = int(row['id'])
        pid = int(row['parent_id'])
        parent_of[nid] = pid
        is_infected_dict[nid] = int(row['is_infected']) == 1
        
        # Запоминаем детей каждого узла для проверки "соседних веток"
        if nid not in children_of:
            children_of[nid] = []
        if pid != -1:
            if pid not in children_of:
                children_of[pid] = []
            children_of[pid].append(nid)
            
    # 2. Восстанавливаем магистральный путь от конца к корню
    path_up = []
    curr = farthest_node_id
    while curr != -1 and curr in parent_of:
        path_up.append(curr)
        curr = parent_of.get(curr, -1)
        
    # 3. Разворачиваем путь, чтобы идти ОТ КОРНЯ ВНИЗ
    path_down = path_up[::-1]
    
    path_nodes = set()
    path_edges = set()
    
    # 4. Идем по стволу и применяем ваши условия
    for i in range(len(path_down)):
        curr = path_down[i]
        
        # УСЛОВИЕ 1: Если сам узел на пути заражен - обрываем ствол
        if is_infected_dict[curr]:
            break
            
        # Добавляем узел в ствол (он точно здоров)
        path_nodes.add(curr)
        
        # Добавляем ребро от родителя
        if i > 0:
            parent = path_down[i-1]
            path_edges.add((curr, parent))
            
        # УСЛОВИЕ 2: Проверяем все исходящие (соседние) ветки
        children = children_of.get(curr, [])
        if len(children) > 0:
            # Если ВСЕ ответвления от этого узла ведут в инфекцию
            all_children_infected = all(is_infected_dict[child] for child in children)
            if all_children_infected:
                break # Обрываем ствол, дальше начинается только красная зона
                
    return path_nodes, path_edges

def plot_trees(data_dir, output_dir, params_str):
    snapshot_files = glob.glob(os.path.join(data_dir, 'snapshot_*.csv'))
    
    def extract_step(filepath):
        match = re.search(r'snapshot_(\d+)\.csv', os.path.basename(filepath))
        return int(match.group(1)) if match else -1

    snapshot_files.sort(key=extract_step)
    
    if not snapshot_files:
        print("[Warning] Snapshot files not found.")
        return

    files_to_plot = snapshot_files[:4]
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14)) 
    
    title = 'Strain Tree Snapshots'
    if params_str: title += f'\n{params_str}'
    fig.suptitle(title, fontsize=16)
    
    axes = axes.flatten()

    for idx, file in enumerate(files_to_plot):
        ax = axes[idx]
        step = extract_step(file)
        df = pd.read_csv(file)

        node_y_dict = calculate_concentrated_layout(df)
        trunk_node_set, trunk_edge_set = find_main_trunk(df)

        coords = {}
        for _, row in df.iterrows():
            nid = int(row['id'])
            coords[nid] = (row['creation_step'], node_y_dict[nid])

        # 3 группы веток (ребер)
        extinct_edges = []
        infected_edges = []
        trunk_edges = []
        
        for _, row in df.iterrows():
            nid = int(row['id'])
            parent_id = int(row['parent_id'])
            is_node_infected = int(row['is_infected']) == 1
            
            if parent_id != -1 and parent_id in coords:
                child_coord = coords[nid]
                parent_coord = coords[parent_id]
                
                v_segment = [parent_coord, (parent_coord[0], child_coord[1])]
                h_segment = [(parent_coord[0], child_coord[1]), child_coord]
                
                # Сортируем ребра по категориям
                if (nid, parent_id) in trunk_edge_set:
                    trunk_edges.append(v_segment)
                    trunk_edges.append(h_segment)
                elif is_node_infected:
                    infected_edges.append(v_segment)
                    infected_edges.append(h_segment)
                else:
                    extinct_edges.append(v_segment)
                    extinct_edges.append(h_segment)

        # 1. Вымершие ветки
        lc_extinct = LineCollection(extinct_edges, colors='blue', alpha=0.3, linewidths=0.5)
        ax.add_collection(lc_extinct)

        # 2. Зараженные ветки
        lc_infected = LineCollection(infected_edges, colors='red', alpha=0.8, linewidths=0.5)
        ax.add_collection(lc_infected)

        # 3. Ствол
        lc_trunk = LineCollection(trunk_edges, colors='green', alpha=1, linewidths=0.7)
        ax.add_collection(lc_trunk)

        df['y_coord'] = df['id'].map(node_y_dict)
        
        is_infected = (df['is_infected'] == 1) & (~df['id'].isin(trunk_node_set))
        is_healthy = (df['is_infected'] == 0) & (~df['id'].isin(trunk_node_set))
        is_trunk = df['id'].isin(trunk_node_set)

        ax.scatter(df[is_healthy]['creation_step'], df[is_healthy]['y_coord'], color='blue', s=2, alpha=0.3) 
        ax.scatter(df[is_infected]['creation_step'], df[is_infected]['y_coord'], color='red', s=2, alpha=1.0) 
        ax.scatter(df[is_trunk]['creation_step'], df[is_trunk]['y_coord'], color='green', s=2, alpha=1.0)

        ax.set_title(f'Step: {step}')
        ax.set_xlabel('Time')
        
        if idx % 2 == 0:
            ax.set_ylabel('Topological Shift (Active -> Up)')
        else:
            ax.set_ylabel('')
            
        if not df.empty:
            ax.set_xlim(-50, df['creation_step'].max() + 50)
            ax.set_ylim(-1, max(node_y_dict.values()) * 1.05 + 1)
            
        if idx == 0:
            # Обновленная легенда (линия + точка)
            legend_elements = [
                Line2D([0], [0], color='#e74c3c', lw=1.5, marker='o', label='Active Branches', markerfacecolor='#e74c3c', markersize=5),
                Line2D([0], [0], color='black', alpha=0.3, lw=1.0, marker='o', label='Extinct Branches', markerfacecolor='black', markersize=4),
                Line2D([0], [0], color='#27ae60', lw=2.0, marker='o', label='Historical Trunk', markerfacecolor='#2ecc71', markersize=5)
            ]
            ax.legend(handles=legend_elements, loc='upper left', fontsize='small')

    for i in range(len(files_to_plot), 4):
        fig.delaxes(axes[i])

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_file = os.path.join(output_dir, 'tree_snapshots.png')
    plt.savefig(out_file, dpi=300, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"Final tree plot saved: {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Visualization of strain tree simulation results.")
    parser.add_argument('-d', '--data_dir', type=str, required=True, help="Path to data directory")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="Path to output directory")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print("Starting plot generation...")
    
    params_str = get_params_string(args.data_dir)
    
    plot_timeseries(args.data_dir, args.output_dir, params_str)
    plot_wave(args.data_dir, args.output_dir, params_str)
    plot_trees(args.data_dir, args.output_dir, params_str)
    
    print("Done!")

if __name__ == "__main__":
    main()
