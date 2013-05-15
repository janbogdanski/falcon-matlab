function w = belka(Q, x)
E = 205000*10^6; %Pa
a = 0.005; %m - bok kwadratowej belki
I = a^4/12; %m4

w  = Q*x^3/(3*E*I);
end

