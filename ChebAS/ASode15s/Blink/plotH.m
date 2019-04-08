for i=1:size(U,1)
   plotEye(0.2,H,i,U,t);
   ylim([-1,1]);
   xlabel(sprintf('t = %f',t(i)));
   
   pause(0.0001);
end