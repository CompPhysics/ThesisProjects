function scaletext(scalefactor)
%SCALETEXT Scales text on the current figure.
%   SCALETEXT(SCALEFACTOR) multiplies the font size of all text
%   in the current figure by SCALEFACTOR.

%   written by Douglas M. Schwarz
%   schwarz@kodak.com
%   5 May 1998

fig = gcf;
sh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')
ax = findobj(fig,'Type','axes');
t = findobj(fig,'Type','text');
for i = 1:length(ax)
    fs = get(ax(i),'FontSize');
    set(ax(i),'FontSize',fs*scalefactor)
end
for i = 1:length(t)
    fs = get(t(i),'FontSize');
    set(t(i),'FontSize',fs*scalefactor)
end
set(0,'ShowHiddenHandles',sh)