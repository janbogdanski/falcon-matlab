
clear;close all;clc
t=linspace(0,0,1000);
count = 1;
figure
tic
while toc < 15
          y=toc;
          plot(t,y,'*');
          
          M(count)=getframe;
          count=count+1;
end
 
movie(M,2,10);