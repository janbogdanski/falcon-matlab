clear;close all;clc

l = 0.1; %m
E = 2e6; %MPa

q = 100; %N/m
a = 0.1; %m - bok kwadratowej belki
I = a^4/12; %m4

m = 0.1; %masa el ruchomego kg

dx = 0.0:.001:l;


h = haptikdevice;

pos = read_position(h);
yold = pos(2) /1000;
ynew = pos(2) /1000;%poczatkowe polozenia, do liczenia pochodnych
told = 0;
count = 1;
figure
tic
while toc < 25
    
    t = toc;
    pos = read_position(h);
    x = pos(1) /1000;
    y = pos(2) /1000;
    ynew = y;
    
    %F = m * ((y - yold)/(toc - told) - yold)/(toc - told)
    dy = y - yold
    dt = toc - told
    v = div(dy,dt)
    a = div(v, dt)
    F = m* a
    %F = 100;
w = 5*F*l^2/(24*E*I) * dx.^2 -  F*l/(12*E*I) * dx.^3

%Q - wartosc sily, obliczona na podstawie m*d2x/dt2 z falcona
%b - punkt uderzenia w belke (belka jeset w srodku, roznica polozen )
%x - zmienna od 0 - l

%Rb = F*x/l;
%Ra = F - Rb;

%Mg = Ra*x - F*(x-b);

          plot(dx,w);
          hold on
          plot(x,y,'*');
          %axis([-100,100,-100,100]) 
          axis([-0.015,0.15,-0.15,0.15]) 
          M(count)=getframe;
          count=count+1;
          yold = ynew;
          told = toc;
end;


close(h);
clear h

movie(M,2,120);