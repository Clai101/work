import os

directory = '.'

old_substring = '-ql'

new_substring = '--ql'

for root, dirs, files in os.walk(directory):
    for filename in files:

        with open(filename, 'r') as file:
            content = file.read()

        new_content = content.replace(old_substring, new_substring)

        with open(filename, 'w') as file:
            file.write(new_content)

print("Замена завершена!")