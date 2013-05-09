function close(h)

if strcmp( get(h.t,'Type') , 'timer')
	stop(h.t);
end

haptik_matlab(6,h.id);


