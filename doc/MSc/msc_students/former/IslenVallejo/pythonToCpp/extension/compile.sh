# compile the interface
echo 'Compiling interface'
swig -c++ -python -I. myModule.i

# compile the sources with 'cpp' extension. 
# Note: -I/... tells the compiler where to 
# find Pithon.h, pyport.h, arrayobject.h
echo 'Compiling sources'
g++ -c -O3 convert.cpp TestCpp.cpp MyArray.cpp myModule_wrap.cxx -I/usr/include/python2.5/

echo 'linking'
# link it together into a shared library so we can use it.
# Note: Because SlaterDeterminant.o includes SlaterMatrix.o,
# it is necesary to link SlaterMatrix.o to the shared library.
g++ -shared -o _myModule.so convert.o myModule_wrap.o 

