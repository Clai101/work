import os
import re
from collections import defaultdict

def read_started_file(file_path):
    """Читает строки из файла started.txt"""
    with open(file_path, 'r') as file:
        return file.readlines()

def generate_file_paths(lines):
    """Генерирует пути к файлам на основе строк из файла started.txt"""
    file_paths = defaultdict(list)
    
    for line in lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split('_')
        if len(parts) < 3:
            print(f"Некорректный формат строки: {line}")
            continue

        type_exp, exp, h = parts[:3]
        d = parts[3] if len(parts) > 3 else ''
        
        file_name = f"{exp}.{h}.{d}.log" if d else f"{exp}.{h}.log"
        directory = os.path.join('.', f"log_{type_exp}")
        file_path = os.path.join(directory, file_name)

        file_paths[type_exp].append(file_path)
        
    return file_paths

def find_files(file_paths):
    """Проверяет наличие файлов и возвращает найденные пути"""
    found_files = defaultdict(list)
    
    for type_exp, paths in file_paths.items():
        for file_path in paths:
            if os.path.isfile(file_path):
                print(f"Файл найден: {file_path}")
                found_files[type_exp].append(file_path)
            else:
                print(f"Файл не найден: {file_path}")
                
    return found_files

def count_basf_statistics(file_path):
    """Подсчитывает количество строк с BASF Execution Statistics и проверяет наличие строки об ошибке чтения файла"""
    basf_count = 0
    skip_file = False

    with open(file_path, 'r') as file:
        for line in file:
            if 'BASF Execution Statistics' in line:
                basf_count += 1
            elif 'This file will not be readable with versions' in line:
                skip_file = True
                break

    return basf_count if not skip_file else 0

def calculate_ratios(found_files):
    """Вычисляет отношение количества BASF Execution Statistics к максимальному значению"""
    max_counts = {'data': 10, 'mc': 100}
    ratios = {}

    for type_exp, files in found_files.items():
        total_basf = sum(count_basf_statistics(file) for file in files)
        max_count = max_counts.get(type_exp, 1)
        ratio = total_basf / max_count
        ratios[type_exp] = ratio
    
    return ratios

def main():
    """Основная функция для выполнения всех шагов"""
    lines = read_started_file('started.txt')
    file_paths = generate_file_paths(lines)
    found_files = find_files(file_paths)
    ratios = calculate_ratios(found_files)
    
    for type_exp, ratio in ratios.items():
        print(f"Тип {type_exp}: отношение фактического количества к максимальному значению = {ratio:.2f}")

if __name__ == "__main__":
    main()
