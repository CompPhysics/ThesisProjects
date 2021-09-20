
from numpy import *

# 0. Tell python where to find the extension module
import sys
sys.path.insert(0, './extension')

# 1. load the wrapped module
import myModule
# 2. define the conversion command
convert = myModule.Convert_MyArray()


# 3. Create some a NxN matrix
N = 5                   # number of rows and columns
A = zeros((N, N))       # zeroes the matrix

# Fill the entries in A
for i in xrange(N):
    for j in xrange(N):
        A[i,j] = i*j

print '\n'
print 'Printing from python'
for i in xrange(N):
    for j in xrange(N):
        print A[i,j]


# Send the numpy 2D-array A (python object) to  C++ 
# 4. Convert the numpy 2D-array object to a C++ MyArray<double> 2D-array
A = convert.py2my(A)

# 5. Call the C++ code
t = myModule.TestCpp(N, A) # Note that N does not need to be converted to a C++ int
t.print_()






