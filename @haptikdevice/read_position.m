function position = read_position(h)
orientation = haptik_matlab(4,h.id);
position = orientation(4,1:3);
