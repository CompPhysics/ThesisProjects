import csv, wla
from carray import *

def read_file(filename, column):
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

def list_to_array(list):
    array = new_carray(len(list))
    for i in range(0,len(list)):   
        carray_set(array,i,list[i]) 
    return array

def fix_length(list):   
    mean = 0.;
    for i in list:
        mean += i;
    mean /= len(list)
                                 #Vector must be of length 2**n. (DWT only)
    n = 1                        #Assume that pasient maintains constant sweating from the time measurements ceased
    lmv = list[len(list)-1]      #until t = 2**n for any n. Add this last measured value at end of vector
    while 2**n < len(list):      
        n += 1
    for x in range(len(list), 2**n):
        list.append(mean)         
    return list, 2**n

def get_arrays(ch_list, v_list, w_list, x_list, stat_list): #Was written for the modwt...(Also works for dwt)
    for i in range(0, len(ch_list)):
        v_list.append(0)
        w_list.append(0)
        x_list.append(0)
    ch = list_to_array(ch_list)
    v = list_to_array(v_list)
    w = list_to_array(w_list)
    x_out = list_to_array(x_list)
    stat = list_to_array(stat_list)
    return ch, v, w, x_out, stat

filename = '/Users/vegard/Master/data/Person'
mean = 0
for i in range(1,31):
    filename2 = filename + str(i)
    if (i == 19):
        filename2 = filename + str(i) + 'c'
    filename3 = 'Person' + str(i)
    channel_number = 4
    j_0 = 2 #2 = Full DWT
#while(channel_number <= 4):
    signal_list = read_file(filename2, channel_number)
    signal_list, array_length = fix_length(signal_list)
    v_list = []
    v_matr = []
    w_list = []
    w_matr = []
    x_list = []
    coeff_list = []
    stat_list = [0,0,0,0,0,0]
    signal, v, w, x_out, stat = get_arrays(signal_list, v_list, w_list, x_list, stat_list)
    object = wla.dwt(4,1)
#array_length = len(signal_list)
    while(array_length >= j_0):#Make different command chains for each algo?
        object.transform(array_length, signal, v, w)
        v_matr.append(v)
        w_matr.append(w)
        array_length /= 2
        coeff_list.append(wla.awc(array_length, w)) #awc should be friend of all classes...
    #Must store old v and w in a list of lists
        for i in range(0,array_length):
            carray_set(signal,i,carray_get(v,i))
#print coeff_list
    coeff = list_to_array(coeff_list)
    wla.linear_regression(len(coeff_list), coeff, stat)
#for i in range (0,6): #r og r_sq tilsynelatende feil selv om den kalkulerer de andre parameterene riktig???
    print filename3, ": ", carray_get(stat, 5)
    mean += carray_get(stat, 5)
mean /= 30
print mean
#channel_number += 1





#signal_list = [1,2,3,4,5,6,7,8,9,10]
#l = 0
#empty_list = []
#signal_list = []
#
#ifile = open('nile.dat','r')
#for value in ifile:
#    signal_list.append(float(value))
#    empty_list.append(0)
#    l+=1
#signal = list_to_array(signal_list)
#v = list_to_array(empty_list)
#w = list_to_array(empty_list)
#w_old = list_to_array(empty_list)
#x_out = list_to_array(empty_list)
#test = wla.modwt(4,1)
#test.transform(l,signal,v,w)
#for i in range(0,l):
#    carray_set(signal,i,carray_get(v,i))
#    carray_set(w_old,i,carray_get(w,i))
#test.transform(l,signal,v,w)
#test.inverse_transform(l,v,w,x_out)
#for i in range(0,l):
#    carray_set(v,i,carray_get(x_out,i))
#test.inverse_transform(l,v,w_old,x_out)
#for i in range(0,10):
#    print carray_get(x_out,i)##Test OK for modwt

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
