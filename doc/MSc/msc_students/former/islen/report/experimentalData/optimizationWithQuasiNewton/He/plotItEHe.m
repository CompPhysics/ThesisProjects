v = load('run1');
it = v(:,1);
a  = v(:,2);
b  = v(:,3);
E  = v(:,4);
plot(it, b);
%legend('\alpha_0 = 1.564,\beta_0 = 0.1342,\alpha_{opt} = 1.8378900818459, \beta_{opt} = 0.3703844012544, E_{min} = -2.8908185137222')
%legend('\alpha_0 = 1.564, \alpha_{opt} = 1.8378900818459')
legend('\beta_0 = 0.1342,\beta_{opt} = 0.3703844012544')