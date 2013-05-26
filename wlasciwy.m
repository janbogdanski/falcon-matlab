%matlabrc;
clear;
close all;clc
eps = 3e-3;


l = 0.05; %m
mesh = 0:0.001:l;

E = 205000*10^6; %Pa

q = 100; %N/m
bok = 0.01; %m - bok kwadratowej belki
I = bok^4/12; %m4
%I = 3.06/100000

m = 0.1; %masa el ruchomego kg

F = 0;
force_max_contact = 1;%N
last_Qf = 0;

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
force_scale = 10000; %tyle razy mniejsza sila generowana niz obliczona

tic
while toc < 25
    
    button = read_button(h);
    if button > 0 
        contact = 0;
        trend = 0;
        trend_in = 0;
        last_trend_in = 0;
    end
        
    
    t = toc;
    pos = read_position(h);
    x = pos(1) /1000;
    y = pos(2) /1000;
    z = pos(3) /1000;

    if (x < l) && (x > 0.03)
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
        
        elseif (abs(y - belka(F,x)) > eps && (yold ~= 0) && y~=0 && (sign(yold) ~= sign(y))&& last_trend_in~=0) 
            %gdy duza predkosc, na razie slabe, bo zostaje w kontakcie po
            %obu stronach( na nowo lapie)
            if contact == 1
                contact = 0;
            else
                contact = 1;
            end
        else
            %null;
        end


        if(contact == 1)

            Q = 3 * y * E*I / (x^3); % y => w, x => l
            Qf = Q/force_scale;
            w = Q*dx.^3/(3*E*I);
            wprim = Q*dx.^2/(2*E*I);
          

            %empirycznie dobrany wspolczynnik, > 1 zapobiega drganiom
            %manipulatora w okolicy pol rownowagi i zapewnia ciaglosc
            %kontaktu z belka
            %zbyt duza wartosc - nie wylapuje odpowienio kontaktu
            if (abs(belka(Q,x)) < 1.044*eps)

                if(last_trend_in ~= 0  && last_trend_in == -trend) || last_trend_in == 0

                    contact = 0;
                    konto = 0;
                    [abs((belka(Q,x))), last_trend_in, trend, (last_trend_in == -trend)];
                end
            else

                if((last_Qf - Qf) > force_max_contact)
                    
                    if abs(Qf) < 5
                        f_step = Qf/10;
                    elseif abs(Qf) < 10
                         f_step = Qf/15;
                    elseif abs(Qf) < 30
                         f_step = Qf/40;
                    else
                         f_step = Qf/100;
                    end

                    %lagodny przebieg sily od 0 do Q co 0.05*1
                    for f = 0:f_step:Qf
                        apply_force(h,-f);
                    end
                    last_Qf = Qf;
                else
                    apply_force(h,-Qf);
                    last_Qf = Qf;
                end
                
            end

        else
            apply_force(h,0);
            elsekontakt = 0;
            dx = 0:0.001:l;
            w = dx *0;
        end




    else
        apply_force(h,0);
        contact = 0;
        last_trend_in = 0;
        dx = 0:0.001:l;
        w = dx *0;
    end

    
          clf
          hold on
          %plot(dx,w);
                  %dorysuj
        





%prostopadly
  

    [A,B,C] = prosta(dx(end-1), w(end-1), dx(end),w(end));
    add = dx(end):0.001:l;
    dx_ext = [dx,add];
    w_ext = [w,(-C/B - A/B*add)];
   % plot(add, (-C/B - A/B*add));


          plot(x,y,'o');
          
          step = 0.001;
          da = 0:step:bok;
          
Y = zeros(length(dx_ext),length(da));
for i=1:length(da)
    Y(:,i) = w_ext(:)  - da(i);
end
Y = Y';
X = meshgrid(dx_ext,da);
[sizey,sizex] = size(X);
if(sizex > 1 && sizey > 1 && length(X) == length(Y))
    color = zeros(sizey,sizex);
    center = (sizey)/2;
    for i=1:sizex
        for j=1:sizey
            if j < center
                z = center - j;
            else
                z = j - center;
            end
            if(i <= length(dx))
                color(j,i) = (sizex-i)*(1)*z;
            else
                color(j,i) = 0;
            end
                
        end
    end

    if w(end)~=0
        
        colormap(jet(ceil(abs(w(end)*1000))));
    else
        colormap(jet(2));
    end

    surf(X,Y,color,'EdgeColor','none')
    

end

          %axis([-100,100,-100,100]) 
          %axis([-0.02,0.08,-0.15,0.15]) troche mala jeszcze
          plot([0.03,0.03],[-1,1],'r'); %czerwona linia ograniczajaca
          plot([0,0],[-0.015,0.005],'LineWidth',4,'Color','k') % 'zamocowanie'
          axis([-0.01,0.07,-0.06,0.06]) 
          axis square;
          grid on
          text(dx(end) + 0.005,w(end), ['w = ',num2str(w(end)), 'm'],'FontSize',10);
          text(dx_ext(end) + 0.005,w_ext(end)-0.005, ['wk = ',num2str(ceil(w_ext(end)*1000)/1000), 'm'],'FontSize',10);
          %M(count)=getframe; %nie potrzebuje nagrywac
          %count=count+1;
          getframe;
          
          yold = ynew;
          told = toc;
end
          %plot(dx,w);
          %hold on
          %plot(x,y,'*');
          %axis([-100,100,-100,100]) 
          %axis([-0.02,0.08,-0.15,0.15]) 
          
close(h);
clear h

%movie(M,2,12);