v = load('run1');
it = v(:,1);
a  = v(:,2);
b  = v(:,3);
E  = v(:,4);

%  plot(it, E);
%  legend('\alpha_0 = 3.85,\beta_0 = 0.08,\alpha_{opt} = 3.9826281304910, \beta_{opt} = 0.1039459074823, E_{min} = -14.5034420665939')

%  plot(it, a);
%  legend('\alpha_0 = 3.85, \alpha_{opt} = 3.9826281304910')

plot(it, b);
legend('\beta_0 = 0.08,\beta_{opt} = 0.1039459074823')