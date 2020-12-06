epitopes = 'C:/sproteinsonecolumn.txt' #Это эпитопы, которые мы хотим проверить на гомологичность с человеческими белками
lisu = "C:/UniprotHuman.txt" #Это база данных белков человека Uniprot
toproy = open(lisu, "r")
conty = toproy.read()
topro = open(epitopes, "r")
cont = topro.read()
part = ''
slices = []
grand = []
tick = 0
for i in range(len(cont)):
    if i == len(cont)-1:
        slices.append(part)
    elif tick == 1:
        if cont[i].isupper():
            part = ''
            part = part + cont[i]
            tick = tick - 1
    elif cont[i].isupper():
        part = part + cont[i]
    elif tick == 0:
        slices.append(part)
        tick = tick + 1
print(len(slices))
for i in range(len(slices)):
    if conty.count(slices[i]) == 0:
        grand.append(slices[i])
    else:
        print(slices[i])
f = open("C:/Tester11.txt", "w")
for i in range(len(grand)):
    f.write(slices[i] + '\n' )
f.close()
print('AhoY!')
 
