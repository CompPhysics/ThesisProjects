function plot_benchmarks(hostname)
% function plot_benchmarks(hostname)

%
% $Id: plot_benchmarks.m,v 2.0 2009/01/22 17:03:38 patrime Exp $
%
% Copyright (c) 2001 Patrick Guio <patrick.guio@fys.uio.no>
%
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%

close all

bench={'daxpy','haney','loop1','loop2','loop3','loop5','loop6','loop8',...
	'loop9','loop10','loop11','loop12','loop13','loop14','loop15','loop16',...
	'loop17','loop18','loop19','loop21','loop22','loop23','loop24','loop25',...
	'loop36','stencil'};

if nargin == 0,
	str=sprintf('blitz-0.9 benchmark on a %s', computer);
else
	str=sprintf('blitz-0.9 benchmark on %s (%s)', hostname, computer);
end
str=sprintf('%s\nCXX=g++ CXXFLAGS= -O3 -funroll-loops -fstrict-aliasing -fomit-frame-pointer -ffast-math', str);
str=sprintf('%s\nF77=ifort FFLAGS= -O3 -Zp16 -ip -pad -unroll -fno-alias -safe_cray_ptr', str);
if length('ifort')
	str=sprintf('%s\nFC=ifort FCFLAGS= -FR -O3 -Zp16 -ip -pad -unroll -fno-alias -safe_cray_ptr', str);
end
str=strrep(str,'_','\_');
h=text(0.5,0.5,str);
set(h,'HorizontalAlignment','center')
set(h,'FontSize',18)
set(h,'FontWeight','demi')
set(gca,'visible','off')
orient landscape
print -dpsc benchmarks.ps

for i=1:length(bench),
	eval(bench{i})
	hs=get(gca,'children')';
	for h=hs, set(h,'linewidth',1.5) , end
	legend
	orient landscape
	print -dpsc -append benchmarks.ps
end


