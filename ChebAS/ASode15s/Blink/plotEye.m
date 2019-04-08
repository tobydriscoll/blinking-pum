function [] = plotEye(pctClosed,H,t_in,U,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

H.sample(U(t_in,1:length(H))');

result = blink(pctClosed,[128 128]);

grid = {result.disc.square.points.x result.disc.square.points.y};

F = H.evalfGrid(grid);

[X,Y]= result.grid(t(t_in));

pcolor(X,Y,F); shading interp;

end

