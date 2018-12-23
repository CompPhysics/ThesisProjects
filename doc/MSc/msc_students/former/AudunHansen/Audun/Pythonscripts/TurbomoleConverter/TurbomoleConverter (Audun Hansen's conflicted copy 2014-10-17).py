#Turbomole converter
#Converts basis files from Basis Set Exchange to class-files readable by Fermion Mingle

filename = "STO6G_H.txt"
f = open(filename)

comments = []

for i in f:
    I = i.split()
    if len(I)==0:
        continue
    if I[0] == '*':
        continue
    if I[0] == '#':        
        comments.append("\\")
        for e in I:
            comments[-1] += e + " "
    #print I
    if I[0] in "123456789":
        print I  
            

print comments

            