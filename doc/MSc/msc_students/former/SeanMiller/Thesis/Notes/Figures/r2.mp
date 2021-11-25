% r2.mp -- do not edit, generated automatically by renorm2.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 10*28.45274pt#, 6*28.45274pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=28.45274pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
straight_galleries;
vleft (__ph1, __i1, __i2);
vright (__ph2, __o1, __o2);
vconnect ("plain", __i1, __v1, __i2);
vconnect ("double,left", __v1, __v2, __v1);
vconnect ("plain", __o1, __v2, __o2);
vdot (__v2);
vdot (__v1);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
