
%vw = VideoWriter('blink.mp4','MPEG-4');
%vw.FrameRate = 30;
%open(vw);
for i=1:size(U,1)
   plotEye(Blinks{1}.percentClosed,H,i,U,t);
   ylim([-1,1]);
   xlabel(sprintf('t = %f',t(i)));
   colorbar();
 %  writeVideo(vw,getframe(gcf));
   pause(0.0001);
end

%close(vw);