text = open('new_AB_model.txt', 'r')

for line in text:
    if 'nan' in line:
        pass
    line  = line.replace('\n', '')
    print line#.replece('\n', '')
