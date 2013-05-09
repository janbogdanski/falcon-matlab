%get a device
clear;close all;clc
h = haptikdevice;

%run simulation for 10 seconds

t=linspace(0,0,100);
count = 1;
figure
tic
while toc < 8
    %read probe position
    pos = read_position(h);
         
      plot(pos(1),pos(2),'*');
          axis([-100,100,-100,100]) 
          M(count)=getframe;
          count=count+1;
        

 


    %check for collision and send back force feedback
     %if pos(2)<0
     %    write(h, -1 * [0 pos(2) 0]);
    % else
    %     write(h,[0 0 0]);
   %  end
    
end

close(h);
clear h

movie(M,2,10);