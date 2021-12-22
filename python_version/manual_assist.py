import os


files_to_open = os.listdir(".")
i = 0
while i + 1 < len(files_to_open):
    j = i + 1
    while j < len(files_to_open):
        if files_to_open[i] == files_to_open[j]:
            files_to_open.pop(j)
            continue
        j+=1

    i+=1


for item in files_to_open:
    if item.count(".m"):
        choice = bool(int(input("Do you want to open " + item + " ?(0 = no, 1 = yes):")))
        if choice:
            os.system("code " + item + " " + item.replace('.m', '.py'))
            option = input("Do you want to continue?(y/n):")
            if option.lower() == 'n':
                break
