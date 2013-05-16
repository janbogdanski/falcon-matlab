matlabrc;clear;close all;clc
eps = 3e-3;


l = 0.05; %m
mesh = 0:0.001:l;

E = 205000*10^6; %Pa

q = 100; %N/m
a = 0.005; %m - bok kwadratowej belki
I = a^4/12; %m4
%I = 3.06/100000

m = 0.1; %masa el ruchomego kg

F = 0;


h = haptikdevice;

pos = read_position(h);
yold = pos(2) /1000;
ynew = pos(2) /1000;%poczatkowe polozenia, do liczenia pochodnych
told = 0;
count = 1;
figure

ret = [];
contact = 0;
trend = 0;
trend_in = 0;
last_trend_in = 0;

tic
while toc < 120
    
    t = toc;
    pos = read_position(h);
    x = pos(1) /1000;
    y = pos(2) /1000;
    z = pos(3) /1000;

    if (x < l) && (x > 0)
        %F = m * ((y - yold)/(toc - told) - yold)/(toc - told)

       dx = 0.0:.001:x; 
       ynew = y;
       
        dy = y - yold;
        dt = toc - told;

        v = div(dy,dt);
        a = div(v, dt);

        F = m * a;

        if(sign(dy) ~=0)
            trend = sign(dy);
        end

        if contact == 0 && sign(dy) ~= 0
            trend_in = sign(dy);
            [trend_in toc];
        end

        %zderzenie
        if(abs(y - belka(F,x)) < eps) && (contact == 0)
            if(last_trend_in ~= 0  && last_trend_in == trend) || last_trend_in == 0

                contact = 1;
                trend_in = trend;
                last_trend_in = trend;

            end
        else
            %contact = 0;
        end


        if(contact == 1)


            Q = 3 * y * E*I / (x^3); % y => w, x => l
            w = Q*dx.^3/(3*E*I);
            [Q, F toc];
            
            


            [sign(y) abs(belka(Q,x)), eps, x,Q];


            if (abs(belka(Q,x)) < eps)

                if(last_trend_in ~= 0  && last_trend_in == -trend) || last_trend_in == 0

                    contact = 0;
                    konto = 0
                    [abs((belka(Q,x))), last_trend_in, trend, (last_trend_in == -trend)];
                end
            else
                apply_force(h,-Q);
            end

            %F = 100;
            %wspornikowa
            %w = 5*F*l^2/(24*E*I) * dx.^2 -  F*l/(12*E*I) * dx.^3;

            %wsp, sila na 'koncu'
            %w = F*dx.^3/(3*E*I);
        else
            apply_force(h,0);
            elsekontakt = 0;
            dx = 0:0.001:l;
            w = dx *0;
        end
%if length(w) >= 1
%    last = w(end-1);
%else
%    last = w(end);
%end
%a =  (-1)*(w(end) - last / 0.1);

%sizeMesh = length(mesh);
%reszta = w(end) + a*(mesh(length(w):length(mesh))-length(w)*0.001);
%w = [w, reszta]

%reszta_dx = mesh(length(dx):length(mesh))
%dx = [dx,reszta_dx]


%dx liczony od zamocowania do miejsca dzialania sily
%nalezy obliczyc kat stycznej i przedluzyc


       % ret = [ret; [dy, dt, v, a, F, 0, x,y,z, w(end)]];

    else
        contact = 0;
        last_trend_in = 0;
        dx = 0:0.001:l;
        w = dx *0;
    end

    


%Q - wartosc sily, obliczona na podstawie m*d2x/dt2 z falcona
%b - punkt uderzenia w belke (belka jeset w srodku, roznica polozen )
%x - zmienna od 0 - l

%Rb = F*x/l;
%Ra = F - Rb;

%Mg = Ra*x - F*(x-b);

          clf
          plot(dx,w);
          hold on
          plot(x,y,'o');
          %axis([-100,100,-100,100]) 
          axis([-0.02,0.08,-0.15,0.15]) 
          M(count)=getframe;
          count=count+1;
          yold = ynew;
          told = toc;
end
          plot(dx,w);
          hold on
          plot(x,y,'*');
          %axis([-100,100,-100,100]) 
          axis([-0.02,0.08,-0.15,0.15]) 
          
close(h);
clear h

%movie(M,2,12);