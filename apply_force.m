function apply_force(handle,F)

F
limit = 5; %max 12!
if(abs(F) < limit)
    force = [0 F, 0];
    write(handle,force);

else
    force = [0 sign(F)*limit, 0];
    write(handle,force)
end