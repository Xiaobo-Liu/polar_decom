function plot_map()
x = -3:0.001:3;
       y = 1/2*x.*(3-x.^2);
       plot(x,y);
       axis([-3 3 -3 3]);
end

