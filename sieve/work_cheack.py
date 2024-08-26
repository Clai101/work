import wandb
import os
import re

def read_started_file(file_path):
    """Читает строки из файла started.txt"""
    with open(file_path, 'r') as file:
        return file.readlines()

def generate_file_paths(lines):
    """Генерирует пути к файлам на основе строк из файла started.txt"""
    file_paths = []
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

        file_paths.append(file_path)
    return file_paths

def find_files(file_paths):
    """Проверяет наличие файлов и возвращает найденные пути"""
    found_files = []
    for file_path in file_paths:
        if os.path.isfile(file_path):
            print(f"Файл найден: {file_path}")
            found_files.append(file_path)
        else:
            print(f"Файл не найден: {file_path}")
    return found_files

def extract_run_number(file_paths):
    """Извлекает значения RunNo и Total Runs из логов и возвращает результаты"""
    pattern_dedx = re.compile(r"INFO\s+:\s+dEdxCalib: ExpNo\s*=\s*\d+\s+RunNo\s*=\s*(\d+)\s+Gain\s*=\s*\d+\.\d+")
    pattern_pntdb = re.compile(r"INFO\s+:\s+PntDB:TPntDB::Get> exp\s*=\s*\d+,\s*run\s*=\s*(\d+),\s*version\s*=\s*1")
    stats_pattern = re.compile(r"^-{2,}\s+BASF Execution Statistics\s+-{2,}")

    results = []

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        found_stats_section = False
        last_match_dedx = None
        last_match_pntdb = None

        for line in lines:
            if stats_pattern.search(line):
                found_stats_section = True

            if found_stats_section:
                match_dedx = pattern_dedx.search(line)
                if match_dedx:
                    last_match_dedx = match_dedx

        last_match_pntdb = None
        for line in lines:
            match_pntdb = pattern_pntdb.search(line)
            if match_pntdb:
                last_match_pntdb = match_pntdb

        run_no_dedx = last_match_dedx.group(1) if last_match_dedx else "0"
        total_runs = last_match_pntdb.group(1) if last_match_pntdb else "Не найдено"

        results.append({
            'file': file_path,
            'run_no_dedx': run_no_dedx,
            'total_runs': total_runs
        })

    return results

def print_results(results):
    """Выводит результаты проверки"""
    print("\nРезультаты:")
    for result in results:
        print(f"Файл: {result['file']}, Run No (dEdxCalib): {result['run_no_dedx']}, Total Runs (PntDB): {result['total_runs']}")

def main():
    """Основная функция для выполнения всех шагов"""
    lines = read_started_file('started.txt')
    file_paths = generate_file_paths(lines)
    found_files = find_files(file_paths)
    results = extract_run_number(found_files)
    print_results(results)

if __name__ == "__main__":
    main()