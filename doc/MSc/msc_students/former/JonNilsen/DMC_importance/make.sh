#!/bin/bash

compiler="g++"
options="-pg -O3"
opt_options="-O3"
lib_options="-pg -O3"
debug="-g"
exec="DMC_importance.app"
cppfiles=`echo $exec | sed 's/\.app/.cpp/g'`
cppfiles="$cppfiles walker.cpp particle.cpp functions.cpp main.cpp"
ofiles=`echo $cppfiles | sed 's/\.cpp/.o/g'`
libos="QickArray/QickArray.o random/random.o"

echo rm -rfv $exec
rm -rfv $exec
echo rm -rfv *.o
rm -rfv *.o

echo cd QickArray
cd QickArray
echo rm -rfv *.o
rm -rfv *.o
echo $compiler $lib_options -c QickArray.cpp
$compiler $lib_options -c QickArray.cpp
echo cd ..
cd ..

echo cd random
cd random
echo rm -rfv *.o
rm -rfv *.o
echo $compiler $lib_options -c random.cpp
$compiler $lib_options -c random.cpp
echo cd ..
cd ..

for file in $cppfiles; do
  echo $compiler $options -c $file
  $compiler $options -c $file
done

echo $compiler $options $libos -o $exec $ofiles
$compiler $options $libos -o $exec $ofiles
