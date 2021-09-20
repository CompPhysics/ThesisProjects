load 'execTimeHeCppNew.data'
cpp = execTimeHeCppNew(:,2);
load 'execTimeHePython.data'
py = execTimeHePython(:,2);
mc = execTimeHeCppNew(:,1);
  diff=abs(cpp-py);

archivo = fopen('relativeTimePyCpp.data', 'w');
  format long e;
for i=1:length(cpp)
	fprintf(archivo, '%d \t%g \t %g \t %g\n',mc(i),cpp(i),py(i),diff(i));

end

fclose(archivo);
