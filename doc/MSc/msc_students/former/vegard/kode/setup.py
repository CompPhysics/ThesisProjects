#!/usr/bin/env python

import commands, os
from distutils.core import setup, Extension

cmd = 'swig -python -c++ carray.i'
print 'Running swig:', cmd
failure, output = commands.getstatusoutput(cmd)
carray_module = Extension('_carray', sources=['carray_wrap.cxx'])
cmd = 'swig -python -c++ wla.i'
#cmd = 'swig -python -c++ testmod.i'
print 'Running swig:', cmd
failure, output = commands.getstatusoutput(cmd)
wla_module = Extension('_wla', sources=['wla_wrap.cxx', 'wla.cpp'])
#testmod_module = Extension('_testmod', sources=['testmod_wrap.cxx', 'testmod.cpp'])

setup (name = 'my_modules',
       version = '0.1',
       author      = "VA",
       #description = """testing""",
       ext_modules = [carray_module, wla_module], 
#       ext_modules = [carray_module, testmod_module]
       )

#python setup.py build_ext --inplace
