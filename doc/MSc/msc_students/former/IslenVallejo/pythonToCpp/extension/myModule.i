%module myModule
%{
#include "MyArray.h"
#include "convert.h"
#include "TestCpp.h"
%}
%include "convert.h"
%include "TestCpp.h"
