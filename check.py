fc = open("CPU")
fg = open("GPU")

for i in fc.readlines():
    for j in fg.readlines():
        if i == j:
           print(i)

fc.close()
fg.close()
