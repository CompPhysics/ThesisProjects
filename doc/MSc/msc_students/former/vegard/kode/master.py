import csv, wla, math, time #time for testing purposes
from carray import *

class Data:
    def __init__(self, option):
        if option is 0:
            self.time = []
            self.etime = []
            self.signal = []
            self.e1 = []
            self.e2 = []
            self.e3 = []
            self.e4 = []
            self.epoch_times = {}
        elif option is 1:
            self.wv1 = []
            self.wv2 = []
            self.wv3 = []
            self.wv4 = []
            self.boundaries = []
        elif option is 2:
            self.ds1 = []
            self.ds2 = []
            self.ds3 = []
            self.ds4 = []
            self.boundaries = []
        elif option is 3:
            self.y1 = []
            self.y2 = []
            self.y3 = []
            self.y4 = []
            self.y_hat1 = []
            self.y_hat2 = []
            self.y_hat3 = []
            self.y_hat4 = []
            self.ci1 = []
            self.ci2 = []
            self.ci3 = []
            self.ci4 = []
            self.out1 = [0,0,0]
            self.out2 = [0,0,0]
            self.out3 = [0,0,0]
            self.out4 = [0,0,0]
        elif option is 4:
            self.hurst1 = []
            self.hurst2 = []
            self.hurst3 = []
            self.hurst4 = []
            self.out1 = [0,0,0]
            self.out2 = [0,0,0]
            self.out3 = [0,0,0]
            self.out4 = [0,0,0]

    def set_time(self, freq):
        self.time = []
        self.etime = []
        for i in range(0, len(self.signal)):
            self.time.append(i/float(freq*60))
        for j in range(0, len(self.e1)):
            self.etime.append(j/float(freq*60))
            
class Settings:
    def __init__(self, option): #Use a parameter as input to set different defaults for different tools?
        if option is 0:
            self.filetype = 1 
            self.channel = 2 
            self.f = 8 
            self.el = 10 
        elif option is 1:
            self.active = 1
            self.transtype = 1
            self.l = 8
            self.type = 2
            self.refl = 0
            self.j_0 = 6
        elif option is 2:
            self.active = 1
            self.transtype = 1
            self.l = 8
            self.type = 2
            self.refl = 0
            self.j_0 = 6
        elif option is 3: 
            self.active = 1
            self.transtype = 1
            self.l = 8
            self.type = 2
            self.refl = 0
            self.j_0 = 6
            self.j_min = 1
            self.j_max = 6
            self.bias = 0
            self.p = 95
        elif(option == 4):
            self.active = 1
            self.transtype = 1
            self.l = 8
            self.type = 2
            self.refl = 0
            self.j_0 = 6
            self.j_min = 1
            self.j_max = 6
            self.bias = 0
            self.smoothlevel = 6

def read_oldfile(filename, column):     #Old filetype
    def comment_stripper(iterator):
        n = 1
        for line in iterator:
            if n > 2:
                if line[:11] == 'Kommentarer':
                    break
                yield line
            n+=1
    def convert_to_double(list):
        return [float(number.split(',')[0]+'.'+number.split(',')[1]) for number in list]
            
    input = csv.reader(comment_stripper(open(filename)), delimiter='\t')
    inputlist = convert_to_double([data[column] for data in input])
    return inputlist

def read_newfile(filename):           #New filetype
    def is_number(x):
        try:
            int(x)
            return True
        except:
            return False
    def convert_to_double(list):
        return [float(number.split(',')[0]+'.'+number.split(',')[1]) for number in list]
    
    input = []
    epoch_times = {}
    time = 1
    tabs = 0

    ifile = open(filename)
    line_list = ifile.readlines()
    ifile.close()
    for datax in line_list[-1]:
        if datax == '\t':
            tabs+=1

    if tabs == 5:
        conductance = 2
        event = 5
    elif tabs == 6:
        conductance = 3
        event = 6
    else:
        print 'File format error'

    for line in line_list:
        data = line.split()
        if len(data) > 0:
            if is_number(data[0]):
                input.append(data[conductance])
                if data[event] != '0':
                    epoch_times[data[event]] = (data[time], len(input)) 
    input = convert_to_double([number for number in input])
    ifile.close()
    return input, epoch_times

def save(ofilename, title, genset, ahurstset, ihurstset, ahurstdata, ihurstdata):
    ofile = open(ofilename, 'a')
    ofile.write(title+"\n")
    if ahurstset.active is 1:
        text = "Averaged estimates: \n"
        ofile.write(text)
        text = "Epoch 1:\t %.5f \n"%(ahurstdata.out1[0])
        ofile.write(text)
        text = "Epoch 2:\t %.5f \n"%(ahurstdata.out2[0])
        ofile.write(text)
        text = "Epoch 3:\t %.5f \n"%(ahurstdata.out3[0])
        ofile.write(text)
        if genset.filetype is 0:
            text = "Epoch 4:\t %.5f \n"%(ahurstdata.out4[0])
            ofile.write(text)
        text = "Variance: \t %.5f \t SD: %.5f \n"%(ahurstdata.out1[1],ahurstdata.out1[2])
        ofile.write(text)
    if ihurstset.active is 1:
        text = "Instant estimates: \n"
        ofile.write(text)
        text = "Epoch 1:\t %.5f \n"%(ihurstdata.out1[0])
        ofile.write(text)
        text = "Epoch 2:\t %.5f \n"%(ihurstdata.out2[0])
        ofile.write(text)
        text = "Epoch 3:\t %.5f \n"%(ihurstdata.out3[0])
        ofile.write(text)
        if genset.filetype is 0:
            text = "Epoch 4:\t %.5f \n"%(ihurstdata.out4[0])
            ofile.write(text)
        text = "Variance:\t %.5f \t SD: %.5f \n"%(ihurstdata.out1[1],ihurstdata.out1[2])
        ofile.write(text)
    ofile.close()

def get_epochtimes(filename, signal, keys, times):
    ifile = open(filename)
    lines = ifile.readlines()
    i = 0
    j = 0
    for line in lines:
        j += 1
        if j > 2:
            time = line.split()[0]
            if(time == times[i]):
                signal.epoch_times[keys[i]] = (time, j)
                i += 1
                if(i > 3):
                    break

    ifile.close()

def get_epochs(settings, signal):   
    freq = settings.f
    time = settings.el*60
    samples = freq*time
    start_one = signal.epoch_times['1'][1]-1
    end_one = start_one + samples
    start_two = signal.epoch_times['2'][1]-1
    end_two = start_two + samples
    start_three = signal.epoch_times['3'][1]-1
    end_three = start_three + samples
    signal.e1 = [signal.signal[i] for i in range(start_one, end_one)]
    signal.e2 = [signal.signal[i] for i in range(start_two, end_two)]
    signal.e3 = [signal.signal[i] for i in range(start_three, end_three)]
    if settings.filetype is 0:
        start_four = signal.epoch_times['4'][1]-1
        end_four = start_four + samples
        signal.e4 = [signal.signal[i] for i in range(start_four, end_four)]
 
def get_maxlevel(transtype, n, l): #Is this bulletproof??
    if transtype is 0:
        l_j = l
        if(n-l_j > 0):
            j = 2
            while l_j <= n:
                l_j = (2**j-1)*(l-1)+1
                j += 1
            j -= 2
            return j
        else:
            return 1
    elif transtype is 1:
        limit = math.log((n/(l-1)+1),2)
        j = 1
        while j < limit:
            j += 1
        j -= 1 #-2????
        return j

def get_ahurstmaxlevel(transtype, n, l, bias):
    if transtype is 0:
        if bias is 0:
            l_j = math.ceil(.5*(l-2))
            n_j = math.floor(.5*n-1)
            if(n-l_j > 0):
                j = 2
                while l_j <= n_j:
                    l_j = math.ceil((l-2)*(1-1/(math.pow(2,j))))
                    n_j = math.floor((n/math.pow(2,j))-1)
                    j += 1
                j -= 2
                return j
            else:
                return 1
        elif bias is 1:
            return math.floor(math.log(n,2))
    if transtype is 1:
        if bias is 0:
            l_j = l
            if(n-l_j > 0):
                j = 2
                while l_j <= n:
                    l_j = (2**j-1)*(l-1)+1
                    print n-l_j
                    j += 1
                j -= 2
                print j
                return j
            else:
                return 1
        elif bias is 1:
            return math.floor(math.log(n,2))
  
def list_to_array(list):
    array = new_carray(len(list))
    for i in range(0,len(list)):   
        carray_set(array,i,list[i]) 
    return array
   
def wa(settings, signal, data):
    object = wla.toolbox()
    transtype = settings.transtype
    l = settings.l
    type = settings.type
    refl = settings.refl
    j_0 = settings.j_0
    n = len(signal.etime) #Assume constant length!!
    print transtype, l, type, refl, j_0, n

    for i in range(0, 2*(j_0+1)):
        data.boundaries.append(0)
    e1 = list_to_array(signal.signal)
    e2 = list_to_array(signal.e2)
    e3 = list_to_array(signal.e3)
    boundaries = list_to_array(data.boundaries)
    wv1 = new_cmatrix(j_0+1,n)
    wv2 = new_cmatrix(j_0+1,n)
    wv3 = new_cmatrix(j_0+1,n)

    object.wa(transtype, l, type, refl, j_0, n, boundaries, e1, wv1)
    object.wa(transtype, l, type, refl, j_0, n, boundaries, e2, wv2)
    object.wa(transtype, l, type, refl, j_0, n, boundaries, e3, wv3)
    #Check for filetype and do extra round...
    if transtype is 0:
        for i in range(j_0+1,0,-1):
            level = i-1
            dummy = 0
            tmp_list1 = []
            for j in range(0,n):
                tmp_list1.append(0)
            for k in range(0,n,2**i):
                tmp_list1[k] = (cmatrix_get(wv1,level,dummy))
                dummy +=1
            data.wv1.append(tmp_list1)
        #boundaries too...
    elif transtype is 1:
        for i in range(j_0+1,0,-1): #Now it plots only the first half of V_j0 
            j = i-1
            tmp_list1 = []
            tmp_list2 = []
            tmp_list3 = []
            for k in range(0,n):
                tmp_list1.append(cmatrix_get(wv1,j,k))
                tmp_list2.append(cmatrix_get(wv2,j,k))
                tmp_list3.append(cmatrix_get(wv3,j,k))
            data.wv1.append(tmp_list1)
            data.wv2.append(tmp_list2)
            data.wv3.append(tmp_list3)

    data.boundaries[0] = carray_get(boundaries, 2*j_0)
    data.boundaries[j_0+1] = carray_get(boundaries, 2*j_0+1)
    j = 1
    for i in range(j_0,0,-1):
        data.boundaries[j] = carray_get(boundaries, i-1)
        j += 1
    j = j_0+2
    for i in range(2*j_0,j_0,-1):
        data.boundaries[j] = carray_get(boundaries, i-1)
        j += 1

def mra(settings, signal, data):
    object = wla.toolbox()
    transtype = settings.transtype
    l = settings.l
    type = settings.type
    refl = settings.refl
    j_0 = settings.j_0
    n = len(signal.e1)
    for i in range(0, 2*(j_0+1)):
        data.boundaries.append(0)
    e1 = list_to_array(signal.e1)
    e2 = list_to_array(signal.e2)
    e3 = list_to_array(signal.e3)
    boundaries = list_to_array(data.boundaries)
    ds1 = new_cmatrix(j_0+1,n)
    ds2 = new_cmatrix(j_0+1,n)
    ds3 = new_cmatrix(j_0+1,n)

    object.mra(transtype, l, type, refl, j_0, n, boundaries, e1, ds1)
    object.mra(transtype, l, type, refl, j_0, n, boundaries, e2, ds2)
    object.mra(transtype, l, type, refl, j_0, n, boundaries, e3, ds3)

    for i in range(j_0+1,0,-1):
        j = i-1
        tmp_list1 = []
        tmp_list2 = []
        tmp_list3 = []
        for k in range(0,n):
            tmp_list1.append(cmatrix_get(ds1,j,k))
            tmp_list2.append(cmatrix_get(ds2,j,k))
            tmp_list3.append(cmatrix_get(ds3,j,k))
        data.ds1.append(tmp_list1)
        data.ds2.append(tmp_list2)
        data.ds3.append(tmp_list3)

    data.boundaries[0] = carray_get(boundaries, 2*j_0)
    data.boundaries[j_0+1] = carray_get(boundaries, 2*j_0+1)
    j = 1
    for i in range(j_0,0,-1):
        data.boundaries[j] = carray_get(boundaries, i-1)
        j += 1
    j = j_0+2
    for i in range(2*j_0,j_0,-1):
        data.boundaries[j] = carray_get(boundaries, i-1)
        j += 1

def ahurst(settings, signal, data):
    object = wla.toolbox()
    transtype = settings.transtype
    l = settings.l
    type = settings.type
    refl = settings.refl
    j_0 = settings.j_0
    j_min = settings.j_min
    j_max = settings.j_max#Must set j_0 wvar needs this
    bias = settings.bias
    p = settings.p
    n = len(signal.e1) #Assume same length of all epochs!!!!
    print l, type, refl, j_0, j_min, j_max, bias, p, n
    
    for i in range(0,j_0):
        data.y1.append(0)
        data.y2.append(0)
        data.y3.append(0)
    for j in range (0, j_0):
        data.y_hat1.append(0)
        data.y_hat2.append(0)
        data.y_hat3.append(0)
    for k in range(0,2*j_0):
        data.ci1.append(0)
        data.ci2.append(0)
        data.ci3.append(0)
    #e4 = list_to_array(data.e4)
    e1 = list_to_array(signal.e1)
    e2 = list_to_array(signal.e2)
    e3 = list_to_array(signal.e3)
    y1 = list_to_array(data.y1)
    y2 = list_to_array(data.y2)
    y3 = list_to_array(data.y3)
    #y4 = list_to_array(data.y4)
    y_hat1 = list_to_array(data.y_hat1)
    y_hat2 = list_to_array(data.y_hat2)
    y_hat3 = list_to_array(data.y_hat3)
    #wvar4 = list_to_array(data.y_hat4)
    ci1 = list_to_array(data.ci1)
    ci2 = list_to_array(data.ci2)
    ci3 = list_to_array(data.ci3)
    #ci4 = list_to_array(data.ci4)
    out1 = list_to_array(data.out1)
    out2 = list_to_array(data.out2)
    out3 = list_to_array(data.out3)
    #out4 = list_to_array(data.out4)

    #if settings.filetype is 1:
    #object.ahurst(transtype, l, type, bias, refl, n, j_0, j_min, j_max, p, e1, y1, y_hat1, ci1, out1)#What about j_0???
    #object.ahurst(transtype, l, type, bias, refl, n, j_0, j_min, j_max, p, e2, y2, y_hat2, ci2, out2)
    #object.ahurst(transtype, l, type, bias, refl, n, j_0, j_min, j_max, p, e3, y3, y_hat3, ci3, out3)

    #if settings.filetype is 0:
    #    object.ahurst(transtype, l, type, bias, refl, n, j_min, j_max, p, e4, y4, wvar4, out4)#Error with wvar

    
    for i in range(0,3):
        data.out1[i] = carray_get(out1,i)
        data.out2[i] = carray_get(out2,i)
        data.out3[i] = carray_get(out3,i)
    #for i in range(0,j_0):
    #    data.y1[i] = carray_get(y1,i)
    #    data.y2[i] = carray_get(y2,i)
    #    data.y3[i] = carray_get(y3,i)
    #for i in range(j_min-1, j_max):
    #    data.y_hat1[i] = carray_get(y_hat1,i)
    #    data.y_hat2[i] = carray_get(y_hat2,i)
    #    data.y_hat3[i] = carray_get(y_hat3,i)
    #del data.y_hat1[j_max:]
    #del data.y_hat2[j_max:]
    #del data.y_hat3[j_max:]
    #del data.y_hat1[0:j_min-1]
    #del data.y_hat2[0:j_min-1]
    #del data.y_hat3[0:j_min-1]
    #for i in range(0,2*j_0):
    #    data.ci1[i] = carray_get(ci1,i)
    #    data.ci2[i] = carray_get(ci2,i)
    #    data.ci3[i] = carray_get(ci3,i)

def ihurst(settings, signal, data):
    object = wla.toolbox()
    l = settings.l
    type = settings.type
    refl = settings.refl 
    j_0 = settings.j_0
    j_min = settings.j_min
    j_max = settings.j_max
    n = len(signal.e1)
    print l, type, refl, j_0, j_min, j_max, n
    for i in range(0, n):
        data.hurst1.append(0)
        data.hurst2.append(0)
        data.hurst3.append(0)

    e1 = list_to_array(signal.e1)
    e2 = list_to_array(signal.e2)
    e3 = list_to_array(signal.e3)
    hurst1 = list_to_array(data.hurst1)
    hurst2 = list_to_array(data.hurst2)
    hurst3 = list_to_array(data.hurst3)
    out1 = list_to_array(data.out1)
    out2 = list_to_array(data.out2)
    out3 = list_to_array(data.out3)
    
    object.ihurst(l, type, refl, j_0, j_min, j_max, n, e1, hurst1, out1)
    object.ihurst(l, type, refl, j_0, j_min, j_max, n, e2, hurst2, out2)
    object.ihurst(l, type, refl, j_0, j_min, j_max, n, e3, hurst3, out3)

    for i in range(0,3):
        data.out1[i] = carray_get(out1,i)
        data.out2[i] = carray_get(out2,i)
        data.out3[i] = carray_get(out3,i)
    for i in range(0,n):
        data.hurst1[i] = carray_get(hurst1,i)
        data.hurst2[i] = carray_get(hurst2,i)
        data.hurst3[i] = carray_get(hurst3,i)

def shurst(data):
    l = 8
    type = 2
    bias = 0
    refl = 0
    j_0 = 11
    j_min = 5
    j_max = 11
    n = len(data.e2)
    n_b = 2
    signal = list_to_array(data.e2)
    object = wla.toolbox()
    object.shurst(l, type, refl, j_0, j_min, j_max, n, n_b, signal)

def analysis():
    object = wla.toolbox()
    signal_list = []#Set signal as "global" here
    y_list = [] #Get these as input
    yhat_list = []#Convert to arrays here
    ci_list = []
    output_list = [0,0,0]
    hurst_list = []
    boundaries_list = []
##
#
#    for value in open('/Users/vegard/Master/Hurst05.dat'): #('/Users/vegard/Master/ocean2.dat'):
#        signal_list.append(float(value.split()[1])) #For reading Ingve's files...
#        hurst_list.append(0)
#       	#signal_list.append(float(value))
#    for i in range(j_min-1,j_max):
#        y_list.append(0)
#    for j in range(0, 3*j_max):
#        wvar_list.append(0)
#    for k in range(0, 2*(j_max+1)):
#        boundaries_list.append(0)
#    n = len(signal_list)
#    signal = list_to_array(signal_list)
#    y = list_to_array(y_list)
#    wvar = list_to_array(wvar_list)
#    output = list_to_array(output_list)
#    hurst = list_to_array(hurst_list)
#    boundaries = list_to_array(boundaries_list)
#    wv = new_cmatrix(j_max+1,n)
#    ds = new_cmatrix(j_max+1,n)
#    object = wla.toolbox()
#    #
#
#    if(wa.get_active() == 1):
#        transtype = wa.get_transtype()
#        l = wa.get_l()
#        type = wa.get_type()
#        refl = wa.get_refl()
#        j_max = wa.get_jmax()
#
#        #Get epochs
#        object.wa(transtype, l, type, refl, j_max, n, boundaries, signal, w)
#    
#    if(mra.get_active() == 1):
#        transtype = mra.get_transtype()
#        l = mra.get_l()
#        type = mra.get_type()
#        refl = mra.get_refl()
#        j_max = mra.get_jmax()
#        
#        #Get epochs
#        object.mra(transtype, l, type, refl, j_max, n, boundaries, signal, ds)
#    
#    if(ahurst.get_active() == 1):
#        transtype = ahurst.get_transtype()
#        l = ahurst.get_l()
#        type = ahurst.get_type()
#        refl = ahurst.get_refl()
#        j_min = ahurst.get_jmin()
#        j_max = ahurst.get_jmax()
#        bias = ahurst.get_bias()
#        p = ahurst.get_p()
#    
#        #Get epochs
#        object.ahurst(transtype, l, type, bias, refl, n, j_min, j_max, p, signal, y, wvar, output) #
#
#    if(ihurst.get_active() == 1):
#        l = ihurst.get_l()
#        type = ihurst.get_type()
#        refl = ihurst.get_refl() #Not implemented yet in C++
#        j_min = ihurst.get_jmin()
#        j_max = ihurst.get_jmax()
#    
#        #Get epochs
#        object.ihurst(l, type, j_min, j_max, len(signal_list), signal, hurst, output)
#        
        
    transtype = 1
    l = 8
    type = 2
    bias = 0
    refl = 0
    j_0 = 6
    j_min = 2
    j_max = 6
    p = 95
    
    for value in open('/Users/vegard/Master/ecg.dat'): #('/Users/vegard/Master/ocean2.dat'):
        #signal_list.append(float(value.split()[1])) #For reading Ingve's files...
        hurst_list.append(0)
      	signal_list.append(float(value))
    for i in range(0,j_0):
        y_list.append(0)
    for j in range(0, j_0):
        yhat_list.append(0)
    for k in range(0,2*j_0):
        ci_list.append(0)
    for k in range(0, 2*(j_max+1)):
        boundaries_list.append(0)
    n = len(signal_list)
    signal = list_to_array(signal_list)
    y = list_to_array(y_list)
    yhat = list_to_array(yhat_list)
    ci = list_to_array(ci_list)
    output = list_to_array(output_list)
    hurst = list_to_array(hurst_list)
    boundaries = list_to_array(boundaries_list)
    wv = new_cmatrix(j_max+1,n)
    #ds = new_cmatrix(j_max+1,n)
    #object = wla.toolbox()
   #For iHurst use j_min:1, j_max:10
    object.ahurst(transtype, l, type, bias, refl, n, j_0, j_min, j_max, p, signal, y, yhat, ci, output) 
    #object.ihurst(l,type,refl,j_min,j_max,n,signal,hurst,output) #Too high j_max will cause segmentation fault
    #object.wa(transtype, l, type, refl, j_max, n, boundaries, signal, wv)
    #print "H =", carray_get(output, 0)
    #print "var:", carray_get(output, 1), " SD Hurst", carray_get(output,2)
    #object.mra(transtype, l, type, refl, j_max, n, boundaries, signal, ds)

if __name__ == "__main__":
    #Declare variables ds, signal etc here... Make set and get funsctions for j_0 etc for the gui.
    #For example get_signal or get_ds should return lists to gui-script for visualization
    print 'This is the wavelet analysis toolbox module. \nRun wlagui.py to use the tools.'
    start = time.time()
    
    signal_list = []
    #for value in open('/Users/vegard/Master/ecg.dat'): #('/Users/vegard/Master/ocean2.dat'):
    #    signal_list.append(float(value.split()[1]))
        #signal_list.append(float(value))
    l = 8
    type = 2
    bias = 0
    refl = 0
    j_0 = 6
    j_min = 1
    j_max = 6
    n = len(signal_list)
    n_b = 6
    signal = list_to_array(signal_list)
    object = wla.toolbox()
#    object.shurst(l, type, refl, j_0, j_min, j_max, n, n_b, signal)

    analysis()
    #ahurstdata = Data()
    #ihurstdata = Data()
    #genset = Settings()
    #ahurstset = Settings()
    #ihurstset = Settings()
    #test = [1,2,3]
    #tittel = "Testfil"
    #ofilename = "/Users/vegard/Master/saves/testfile"
    #ahurstdata.out1 = test
    #ahurstdata.out2 = test
    #ahurstdata.out3 = test
    #save(ofilename, tittel, genset, ahurstset, ihurstset, ahurstdata, ihurstdata)
    #signal_list = []
    #hurst_list = []
    #var_list = []


    #signal_list = read_file('/Users/vegard/Master/data/Person6', 2)
    
    #for value in open('/Users/vegard/Master/Hurst05.dat'): #('/Users/vegard/Master/ocean2.dat'):
    #    signal_list.append(float(value.split()[1]))
    #    hurst_list.append(0)
    #    var_list.append(0)
    #    signal_list.append(float(value.split()[1])) #For reading Ingve's files...
    #signal_list = reflection(signal_list)
  
    #modwt_hurst(8, 2, 1, 17, 1, signal_list)
    #signal_list=reflection(signal_list)
    #signal = list_to_array(signal_list)
    #hurst = list_to_array(hurst_list)
    #var = list_to_array(var_list)
    #l = len(signal_list)
    #object = wla.modwt(12,2)
    #object.iwlse(1,12,l,signal,hurst,var) #12 for Ingves files
    #dwt_mra(4,1,6,signal_list)
    #ofile = open('ht.dat','w')
    #for i in range(0,len(signal_list)/2):
    #    ofile.write(str(carray_get(hurst,i))+'\n')
    
    

    #Length of filter which gives a stabile H is related to the order of the backward difference in the non-stationary
    #signal (A good idea is to add more Daubechies filters to the filter bank)
    #I suspect the parameter d to be important! (A person with abnormal SA would have a much larger d than a normal
    #person)
    #NB!!!!! Check chapter 9.5 against my linear regression model!!!
    #print '-------------------'
    #dwt_hurst(4, 1, 1, 9, 1, sl)

    elapsed = (time.time() - start)
    print 'Time : ', elapsed


        
#Classes:
#GUI class
#Hurst class
#MRA class
#???


############### OLD COMMAND CHAINS ############
#def reflection(list):#Delete
#    copy = list[:]
#    reflected = copy[:]
#    reflected.reverse()
#    copy.extend(reflected)
#    
#    return copy
#
#def slice(min, max, list):#Delete
#    return list[min:max] #min -1???
#
#def dwt_mra(l, family, j_max, signal_list):#Delete
#    object = wla.dwt(l, family)
#    length = len(signal_list)
#    signal = list_to_array(signal_list)
#    ds = new_cmatrix(j_max+1, length)
#    object.mra(length, j_max, signal, ds)
#
#General MRA
#Delete
#def mra(type, family, size, j_0, signal_list): #Reflection or circularity should be an input parameter
#    if(type == 1): #Write two different functions like for Hurst???
#        object = wla.dwt(size, family) #Changed family,size to size,family
#        print 'DWT!!!'
#    elif(type == 2):
#        object = wla.modwt(size, family) #Changed here too
#    
#    length = len(signal_list)
#    signal = new_carray(length)
#    for i in range(0,length):
#        carray_set(signal, i, signal_list[i])
#    ds = new_cmatrix(j_0+1,length)
#    object.mra(length, j_0, signal, ds)
#    
#    ds_list = []
#    for i in range(0,j_0+1):
#        tmp = []
#        for j in range(0, length):
#            tmp.append(cmatrix_get(ds, i, j))
#        ds_list.append(tmp)
#
#    return ds_list
#Delete
#def dwt_hurst(l, family, j_min, j_max, bias, signal_list):
#    object = wla.dwt(l, family)       #Wavevar-function should also return a confidence interval of the estimator
#    v_list = []
#    w_list = []
#    coeff_list = []
#    stat_list = [0,0,0,0,0,0]
#    
#    n = 0            
#    while 2**n < len(signal_list):      
#        n += 1
#    for x in range(len(signal_list), 2**n):
#        signal_list.append(0)                #Pad with zeros
#    signal_length = 2**n
#
#    for i in range(0, signal_length):
#        v_list.append(0)
#        w_list.append(0)
#    
#    signal = list_to_array(signal_list)
#    v = list_to_array(v_list)
#    w = list_to_array(w_list)
#    stat = list_to_array(stat_list)
#
#    for i in range(0, j_max):
#        object.transform(signal_length, signal, v, w)
#        if(i >= j_min):
#            coeff_list.append(object.wavevar(bias, i+1, signal_length, w))
#        for j in range(0, signal_length):
#            carray_set(signal, j, carray_get(v, j))  
#
#    coeff = list_to_array(coeff_list)
#    wla.linear_regression(len(coeff_list), coeff, stat) #Vil helst skrive object.linear_reg...
#    print 'beta_0 : ', carray_get(stat,0)
#    print 'beta_1 : ', carray_get(stat,1)
#    print 'sigma : ', carray_get(stat,2)
#    print 'r : ', carray_get(stat,3)
#    print 'r_sq : ', carray_get(stat,4)
#    print 'H : ', carray_get(stat, 5)
#    
#    #Delete
#def modwt_hurst(l, family, j_min, j_max, bias, signal_list):
#    object = testmod.modwt(l, family)              
#    v_list = []                                
#    w_list = []
#    coeff_list = []
#    conf_list = [0,0]
#    stat_list = [0,0,0,0,0,0]                  
#    signal_length = len(signal_list)
#    for i in range(0, signal_length): #Making v and w too long may cause errors... Need only to be signal_length/2
#        v_list.append(0)
#        w_list.append(0)
#
#    signal = list_to_array(signal_list)
#    v = list_to_array(v_list)
#    w = list_to_array(w_list)
#    conf = list_to_array(conf_list)
#    stat = list_to_array(stat_list)
#
#    for i in range(0, j_max): #Have changed wla.cpp. This syntax makes it (0,1,...,j_max-1). New it (1,2,...,j_max)
#        object.transform(signal_length, i+1, signal, v, w)
#        #for k in range(0,6):
#        #    print carray_get(w, k)
#        #print "-------------"
#        if(i >= j_min):
#            coeff_list.append(object.wavevar(bias, i+1, signal_length, 95, w, conf)) ##Then i, not i+1
#        for j in range(0, signal_length):
#            carray_set(signal, j, carray_get(v, j))  
#
#    coeff = list_to_array(coeff_list)
#    object.wlse(j_min,j_max,signal_length,coeff)
#    #wla.linear_regression(len(coeff_list), coeff, stat) #Vil helst skrive object.linear_reg...
#    #print 'beta_0 : ', carray_get(stat,0)               #I may do this by code modwt::linreg which just passes the
#    #print 'beta_1 : ', carray_get(stat,1)               #parameters to linear_regression. Also make a dwt::linreg 
#    #print 'sigma : ', carray_get(stat,2)                #which do the same thing
#    #print 'r : ', carray_get(stat,3)
#    print 'r_sq : ', carray_get(stat,4)
#    print 'H : ', carray_get(stat, 5)


############### DWT TRANSFORM ############
#    length = len(signal_list)
#    for i in range(0,length/2-10):
#        v_list.append(0)
#        w_list.append(0)
#    object = wla.dwt(4, 1)
#    
#    signal = list_to_array(signal_list)
#    v = list_to_array(v_list)
#    w = list_to_array(w_list)

    #for i in range(0,1):
#    for i in range(0,2):
#        object.transform(length, signal, v, w)
#        for j in range(0, length):
#            carray_set(signal, j, carray_get(v, j))
#        length /= 2
#        print '__________'
#        for k in range(0,15):
#            print carray_get(w,k)


############## SHIFTS #############
#signal_list = []
#    v_list = []
#    w_list = []
#    
#    #signal_list = read_file('/Users/vegard/Master/data/Person6', 2)
#    
#    for value in open('/Users/vegard/Master/ecg.dat'): #('/Users/vegard/Master/ocean2.dat'):
#        signal_list.append(float(value))
#        #signal_list.append(float(value.split()[1])) #For reading Ingve's files...
#    #signal_list = reflection(signal_list)
#  
#    #modwt_hurst(8, 2, 1, 17, 1, signal_list)
#    for i in range(0, len(signal_list)): #Making v and w too long may cause errors...
#        v_list.append(0)
#        w_list.append(0)
#    signal = list_to_array(signal_list)
#    v = list_to_array(v_list)
#    w = list_to_array(w_list)
#    #dwt_mra(4,1,6,signal_list)
#    object = wla.modwt(8,2)
#    object.transform(len(signal_list), 1, signal, v, w)
#    for j in range(0, len(signal_list)):
#        carray_set(signal, j, carray_get(v, j)) 
#    #test = object.h_shift(1)
#    object.shift(0,1,len(signal_list),w)
#    for i in range(0,10):
#        print carray_get(w, i)
#    print '_____'
#    for i in range(len(signal_list)-10, len(signal_list)):
#        print carray_get(w, i)
#    for j in range(0, len(signal_list)):
#        carray_set(signal, j, carray_get(v, j)) 
#    object.transform(len(signal_list), 2, signal, v, w)
#    object.shift(0,2,len(signal_list),w)
#    print '____'
#    print '____'
#    for i in range(0,10):
#        print carray_get(w, i)
#    print '_____'
#    for i in range(len(signal_list)-10, len(signal_list)):
#        print carray_get(w, i)
#    #for i in range(1,6):
#    #    print object.h_shift(i)


#test_hurst_dwt()            
    #test_transform_modwt()
    #test_mra()
    #test()#Running mra command chain in c++ approx 10 times faster than python (using carray_get/set all the time)
    #test_dwt_mra()

#def fix_length(list):   
#    mean = 0.;
#    for i in list:
#        mean += i;
#    mean /= len(list)
#                                 #Vector must be of length 2**n. (DWT only)
#    n = 1                        #Assume that pasient maintains constant sweating from the time measurements ceased
#    lmv = list[len(list)-1]      #until t = 2**n for any n. Add this last measured value at end of vector
#    while 2**n < len(list):      
#        n += 1
#    for x in range(len(list), 2**n):
#        list.append(mean)         
#    return list, 2**n
#
#def get_arrays(signal_list, v_list, w_list, stat_list): #Was written for the modwt...(Also works for dwt)
#    for i in range(0, len(signal_list)): #Not needed?????????????????????
#        v_list.append(0)      #This function is to be deleted
#        w_list.append(0)
#    signal = list_to_array(signal_list)
#    v = list_to_array(v_list)
#    w = list_to_array(w_list)
#    stat = list_to_array(stat_list)
#    return signal, v, w, stat
    

#from copy import deepcopy #Do I need this???


#signal_list = [1,2,3,4,5,6,7,8,9,10]
#l = 0
#empty_list = []
#signal_list = []
#stat_data = list_to_array(stat);
#test = wla.dwt(4,1)
#v_in_list = [1,2,3,4,5,6,7,8]
#v_list = []
#w_list = []
#x_out_list = []
#for i in range(0,array_length/2):
#    v_list.append(0)
#    w_list.append(0)
#for i in range(0,array_length):
#    x_out_list.append(0)
#x_out_list = [0,0,0,0,0,0,0,0]
#v_in = list_to_array(v_in_list)
#v = list_to_array(v_list)
#w = list_to_array(w_list)
#x = list_to_array(v_in_list)
#x_out = list_to_array(x_out_list)
#ecg_list = []
#ifile = open('ecg.dat','r')
#for value in ifile:
#    ecg_list.append(float(value))
#array_length = 2048
#ecg = list_to_array(ecg_list)
#test.transform2(8, v_in, v, w)
#test.transform(8,x)
#for i in range(0,4):
#    carray_set(x,i,carray_get(v,i))
#    carray_set(x,i+4,carray_get(w,i))
#test.inverse_transform(8,x)

#test.inverse_transform2(4096, v, w, x_out)
#test.transform(len(v_in_list),x)
#j = 0
#for i in range(0,4):
#    print carray_get(v,i)
#for i in range(0,4):
#    print carray_get(w,i)
#for i in range(0,8):
#    print carray_get(x,i)
#test2 = wla.modwt(4,1)
#test3 = wla.dwt(4,1)

#ch, array_length = fix_length(ch)

#    test.transform(array_length, ecg, v, w)
    #test.transform(array_length, cpparray)
#    array_length /=2
#    coeff.append(wla.awc(array_length, w))
#    for i in range(0,array_length):
        #carray_set(w,i,0)
#        carray_set(ecg,i,carray_get(v,i))
    #for i in range(0,5):
    #    print carray_get(w,i)
    #print "---"
#for i in range(0,5):
#    print carray_get(v,i)
#print coeff
#ccoeff = list_to_array(coeff)
#wla.linear_regression(len(coeff), ccoeff, stat_data)
#for i in range (0,6): #r og r_sq tilsynelatende feil selv om den kalkulerer de andre parameterene riktig???
#    print carray_get(stat_data, i)
#channel_number += 1

#
#test_* are command chains for testing things...
#

#def test_hurst_dwt(): #Make a c++ function to save all the array conversions... Output: stat_data
    #filename = '/Users/vegard/Master/data/Person'
    #mean = 0

    #for i in range(1,31):
    #    filename2 = filename + str(i)
    #    if (i == 19):
    #        filename2 = filename + str(i) + 'c'
    #    filename3 = 'Person' + str(i)
    #    channel_number = 2
    #    j_0 = 2 #2 = Full DWT
#while(channel_number <= 4):
    #   signal_list = read_file(filename2, channel_number)
    #    signal_list, array_length = fix_length(signal_list)
    #    v_list = []
    #    w_list = []
    #    w_matr = []
    #    x_list = []
    #    coeff_list = []
    #    stat_list = [0,0,0,0,0,0]
    #    signal, v, w, x_out, stat = get_arrays(signal_list, v_list, w_list, x_list, stat_list)
    #    object = wla.dwt(4,1)
#array_length = len(signal_list)
        
        
     #   while(array_length >= j_0):#Make different command chains for each algo?
     #       object.transform(array_length, signal, v, w)
     #       w_matr.append(w)
     #       array_length /= 2
     #       coeff_list.append(wla.awc(array_length, w)) #awc should be friend of all classes...
    #Must store old v and w in a list of lists
     #      for i in range(0,array_length):
     #           carray_set(signal,i,carray_get(v,i))
#print coeff_list
     #   coeff = list_to_array(coeff_list)
     #   wla.linear_regression(len(coeff_list), coeff, stat)
#for i in range (0,6): #r og r_sq tilsynelatende feil selv om den kalkulerer de andre parameterene riktig???
     #   print filename3, ": ", carray_get(stat, 5), carray_get(stat, 3)
     #   mean += carray_get(stat, 5)
    #mean /= 30
    #print 'Mean : ', mean
    #print carray_get(stat,0)
    #print carray_get(stat,1)
    #print carray_get(stat,2)
    #print carray_get(stat,3)
    #print carray_get(stat,4)
    #print coeff_list
    
    #hurstfilename = filename+'31'
    #hurst_list = []
    #for value in open(hurstfilename, 'r'):
    #    hurst_list.append(float(value))
    #hurst_array = list_to_array(hurst_list)
    #array_length = len(hurst_list)
    #test_list=[]
    #v_new_list = []
    #w_new_list = []
    #for i in range(0,array_length):
    #    v_new_list.append(0)
    #    w_new_list.append(0)
    #v_new = list_to_array(v_new_list)
    #w_new = list_to_array(w_new_list)
    #while(array_length >= 2):
        #object.transform(array_length, hurst_array, v_new, w_new)
        #array_length /= 2
        #test_list.append(wla.awc(array_length, w_new))#What is N???
        #for i in range(0,array_length):
        #    carray_set(hurst_array,i,carray_get(v_new,i))      

    #test = list_to_array(test_list)
    #wla.linear_regression(len(test_list), test, stat)
    #print len(signal_list)
    #print carray_get(stat,5), carray_get(stat,3)
    #print test_list

    #####Ecg wavelet coeffisients##########
    #ecg_list = []
    #for value in open('/Users/vegard/Master/ecg.dat'):
    #    ecg_list.append(float(value))
    #array_length = len(ecg_list)
    #ecg = list_to_array(ecg_list)
    #v_new_list = []
    #w_new_list = []
    #for i in range(0,array_length):
    #    v_new_list.append(0)
    #    w_new_list.append(0)
    #v_new = list_to_array(v_new_list)
    #w_new = list_to_array(w_new_list)
    #while(array_length >= 2):
    #    print '-------------------'
    #    object.transform(array_length, ecg, v_new, w_new)
    #    for i in range(0,10):
    #        print carray_get(w_new,i)
    #    array_length /= 2
    #    for i in range(0,array_length):
    #        carray_set(ecg,i,carray_get(v_new,i))  

#    object2 = wla.modwt(4,1)
#    ocean_list = []
#    c_list = []
#    stat_list = [0,0,0,0,0,0]
#    for value in open('/Users/vegard/Master/data/Person31'): #('/Users/vegard/Master/ocean2.dat'): 
#        ocean_list.append(float(value))
#    array_length = len(ocean_list)
#    ocean = list_to_array(ocean_list)
#    v_new_list = []
#    w_new_list = []
#    for i in range(0,array_length):
#        v_new_list.append(0)
#        w_new_list.append(0)
#    v_new = list_to_array(v_new_list)
#    w_new = list_to_array(w_new_list)
#    stat = list_to_array(stat_list)
#    
#    for i in range(0,4):
#        object2.transform(array_length, i, ocean, v_new, w_new)
#        c_list.append(object2.wavevar(0,i+1,array_length,w_new))
#        #print object2.wavevar(0,i+1,array_length,w_new)
#        for k in range(0,array_length):
#            carray_set(ocean,k,carray_get(v_new,k))  
#
#    c = list_to_array(c_list)
#    wla.linear_regression(len(c_list), c, stat) #Vil helst skrive object2.linear_reg...
#    print 'beta_0 : ', carray_get(stat,0)
#    print 'beta_1 : ', carray_get(stat,1)
#    print 'sigma : ', carray_get(stat,2)
#    print 'r : ', carray_get(stat,3)
#    print 'r_sq : ', carray_get(stat,4)
#    print 'H : ', carray_get(stat, 5)
    
        

#def test_transform_modwt():
    #What about calculating H for parts of the signals (for example running only)
#    j_0 = 4
#    channel_number = 2
#    filename = '/Users/vegard/Master/data/Person7'
#    signal_list = []
#    v_list = []
#    w_list = []
#    w_matr = []
#    x_list = []
#    stat_list = [0,0,0,0,0,0]

    #ifile = open('nile.dat','r')
    #for value in ifile:
    #    signal_list.append(float(value))
    #signal_list = read_file(filename, channel_number)

#    signal, v, w, x_out, stat = get_arrays(signal_list, v_list, w_list, x_list, stat_list)
#    object = wla.modwt(2,1)
#    length = len(signal_list)
#    
#    for i in range(0,j_0):
#        object.transform(length, i, signal, v, w)
#        for j in range(0, length):
#            carray_set(signal, j, carray_get(v,j))
#            w_list[j] = carray_get(w,j)
#        w_matr.append(deepcopy(w_list))
        #print wla.awc(length, w)
#    for i in range(0, len(w_matr)):
#        for j in range(0, 10):
#            print w_matr[i][j]
#        print '---'
#    
#    while(j_0 > 0):
#        j_0 -= 1
#        object.inverse_transform(length, j_0,v, w, x_out)
#        for i in range(0, length):
#            carray_set(v, i, carray_get(x_out,i))
#            carray_set(w, i, w_matr[j_0-1][i])
#
#    for i in range(0, 10):
#        print carray_get(x_out, i)
# 
#def test_mra(): #MRA for modwt Make this a c++ function???
#    j_0 = 4
    #channel_number = 2
    #filename = '/Users/vegard/Master/data/Person18'
#    signal_list = []
#    d_list = []
#    d_matr = []
#    s_list = []
#    v_list = []
#    w_list = []
#    w_matr = []
#    x_list = []
#    stat_list = [0,0,0,0,0,0]
#
#    ifile = open('nile.dat','r')
#    for value in ifile:
#        signal_list.append(float(value))
#        d_list.append(0)
#        s_list.append(0)
    #signal_list = read_file(filename, channel_number)
    #for i in range(0, len(signal_list)):
    #    d_list.append(0)
    #    s_list.append(0)

#    signal, v, w, x_out, stat = get_arrays(signal_list, v_list, w_list, x_list, stat_list)
#    object = wla.modwt(2,1)
#    length = len(signal_list)
#    
#    for i in range(0,j_0):
#        object.transform(length, i, signal, v, w)
#        for j in range(0, length):
#            carray_set(signal, j, carray_get(v,j))
#            w_list[j] = carray_get(w,j)
#        w_matr.append(deepcopy(w_list))
#    
#    for i in range(0, length):
#        carray_set(signal,i,0)#Set signal to zero
#        
#    #Compute details D_j:
#    for i in range(0, j_0):
#        for j in range(0, length):
#            carray_set(w, j, w_matr[i][j])
#        object.inverse_transform(length, i, signal, w, x_out)#Signal should be a list named zero...
#        k = i-1
#        while(k >= 0):
#            for l in range(0, length):
#                carray_set(w, l, carray_get(x_out, l))
#            object.inverse_transform(length, k, w, signal, x_out)
#            k -= 1
#        for m in range(0, length):
#            d_list[m] = carray_get(x_out, m)
#        d_matr.append(deepcopy(d_list))
    #d_matr.reverse()
    #for i in range (0,len(d_matr)):
    #    for j in range(0,10):
    #        print d_matr[i][j]
    #    print '---'
    #Must fix the transform algo (it)
    #Compute Smooth S:
    
    #Do the same thing one more time...(V_j as initial input (179))
#    object.inverse_transform(length, j_0-1, v, signal, x_out)#Signal should be a list named zero...
#    k = j_0-2
#    while(k >= 0):
#        for l in range(0, length):
#            carray_set(v, l, carray_get(x_out, l))
#        object.inverse_transform(length, k, v, signal, x_out)
#        k -= 1
#    for m in range(0, length):
#        s_list[m] = carray_get(x_out, m)
    #d_matr.append(deepcopy(d_list))
    #d_matr.reverse()
#
#    for i in range (0,len(d_matr)):
#        for j in range(0,10):
#            print d_matr[i][j]
#        print '---'
#    for i in range(0, length):
#        signal_list[i] = s_list[i]
    #Build the signal again...
#    for i in range(0, len(d_matr)):
#        for j in range(0, length):
#            signal_list[j] += d_matr[i][j]
#    for i in range(0,10):
#        print signal_list[i]
            
#def test():
#    j_0 = 4
#    signal_list = []
#    ifile = open('nile.dat','r')
#    for value in ifile:
#        signal_list.append(float(value))
#    length = len(signal_list)
#    #signal = list_to_array(signal_list)
#    signal = new_carray(length)
#    for i in range(0,len(signal_list)):
#        carray_set(signal, i, signal_list[i])
#    ds = new_cmatrix(j_0+1,length) 
#    object = wla.modwt(2,1)
#    object.mra(length, j_0, signal, ds)
    #del signal#Does this work???
    #del ds
    #What about not deleting arrays/matrices???

#def test_dwt_mra():
#    j_0 = 6
#    signal_list = []
#    for value in open('ecg.dat','r'):
#        signal_list.append(float(value))
#    length = len(signal_list)
#    signal = new_carray(length)
#    for i in range(0,length):
#        carray_set(signal, i, signal_list[i])
#    ds = new_cmatrix(j_0+1,length)
#    object = wla.dwt(6,4)
#    object.mra(length, j_0, signal, ds)
#    del signal #Return Nul pointer instead?
#    del ds

#def get_x(list):
#    delta_t = 1
#    x = []
#    max = list[0]
#    for i in range(0, len(list)):
#        x.append(i*delta_t)
#        if i > max:
#            max = i
#
#    return x, max

