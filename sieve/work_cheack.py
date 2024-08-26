import os
import time
import wandb
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

def calculate_file_ratios(found_files):
    """Вычисляет отношение количества BASF Execution Statistics к максимальному значению для каждого файла"""
    max_counts = {'data': 10, 'mc': 100}
    file_ratios = {}
    counts_by_type = defaultdict(int)

    for type_exp, files in found_files.items():
        max_count = max_counts.get(type_exp, 1)
        for file in files:
            basf_count = count_basf_statistics(file)
            ratio = basf_count / max_count
            file_ratios[file] = ratio
            counts_by_type[type_exp] += basf_count
    
    return file_ratios, counts_by_type

def calculate_type_percentage(lines, counts_by_type):
    """Вычисляет процентное отношение количества файлов к максимальному значению по типу"""
    max_counts = {'data': 10, 'mc': 100}
    type_counts = defaultdict(int)
    
    for line in lines:
        type_exp = line.split('_')[0]
        type_counts[type_exp] += 1

    type_percentages = {}
    for type_exp, count in type_counts.items():
        max_count = max_counts.get(type_exp, 1) * count
        sum_count = counts_by_type.get(type_exp, 0)
        type_percentages[type_exp] = (sum_count / max_count) * 100
    
    return type_percentages

def log_metrics():
    """Функция для выполнения всех шагов и логирования метрик"""
    lines = read_started_file('started.txt')
    file_paths = generate_file_paths(lines)
    found_files = find_files(file_paths)
    file_ratios, counts_by_type = calculate_file_ratios(found_files)
    
    # Логирование отношения для каждого файла
    for file_path, ratio in file_ratios.items():
        wandb.log({f"ratio/{file_path}": ratio})

    # Логирование процента по типам
    type_percentages = calculate_type_percentage(lines, counts_by_type)
    for type_exp, percentage in type_percentages.items():
        wandb.log({f"percentage/{type_exp}": percentage})

def main():
    """Основная функция с циклом для выполнения задания каждую минуту"""
    # Инициализация W&B
    wandb.init(project="sieve", entity="clai101-hse-university")
    
    try:
        while True:
            log_metrics()
            print("Завершена итерация логирования. Ожидание 1 минуту...")
            time.sleep(60)  # Задержка в 1 минуту
    except KeyboardInterrupt:
        print("Логирование остановлено вручную.")
    finally:
        wandb.finish()  # Завершение сессии W&B

if __name__ == "__main__":
    main()
