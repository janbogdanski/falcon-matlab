function write(h,ft)

if nargin<2
 ft = [0 0 0
       0 0 0];
end

[m,n] = size(ft);

if (m == 1 && n ==6)
	ft = [ft(1:3) 
	      ft(4:6)];
elseif (m == 1 && n ==3) 
	ft = [ft 
	      zeros(1,3)];
elseif (m == 3 && n ==2) 
	ft = ft';
elseif (m == 3 && n ==1) 
	ft = [ft zeros(3,1)]';
end

haptik_matlab(5,ft,h.id);