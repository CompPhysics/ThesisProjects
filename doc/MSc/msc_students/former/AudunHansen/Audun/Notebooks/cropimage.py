from matplotlib.pyplot import *


from os import listdir
from os.path import isfile, join


def cropmage(img):
    #autotrim image size
    if len(img[0,0,:]) == 4 or len(img[0,0,:]) == 3:
        im = img[:,:,0]
        Nx = len(im)-1
        Ny = len(im[0])-1
        nx = 0
        ny = 0
    
        #trim image
        while im[nx].min() == 1:
            nx += 1
        while im[:,ny].min() == 1:
            ny += 1
        
        while im[Nx].min() == 1:
            Nx -= 1
        while im[:,Ny].min() == 1:
            Ny -= 1
            
        nx -= 2
        ny -= 2
        Nx += 2
        Ny += 2
    else:
        Nx = len(img)-1
        Ny = len(img[0])-1
        nx = 0
        ny = 0

        im = img
        
    
        #trim image
        while im[nx].min() == 1:
            nx += 1
        while im[:,ny].min() == 1:
            ny += 1
        
        while im[Nx].min() == 1:
            Nx -= 1
        while im[:,Ny].min() == 1:
            Ny -= 1
            
        nx -= 2
        ny -= 2
        Nx += 2
        Ny += 2


    return img[nx:Nx, ny:Ny, :]

def trim_all():
    onlyfiles = [ f for f in listdir("/Users/kinealicegulbrandsen/Dropbox/master thesis audun skau hansen/Notebooks") if isfile(join("/Users/kinealicegulbrandsen/Dropbox/master thesis audun skau hansen/Notebooks",f)) ]
    for f in onlyfiles:
        exte = f[-3:]
        if f[:4] == "CCDT" and exte == "png":
            print f
            img = imread(f)
            img = cropmage(img)
            imsave(f, img)
        


trim_all()

#img = imread("CCSDT_t3_9_29_0.png")
#print len(img[0,0])
#img = cropmage(img)
#imshow(img)
#show()
#img = cropmage(img)
#imsave("CCSDT_t3_9_29_0.png", img)
