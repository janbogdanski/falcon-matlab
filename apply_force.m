function apply_force(handle,F)
F;
if(abs(F) < 4)
    force = [0 F, 0];
   % write(handle,force);

else
    force = [0 sign(F)*2, 0];
    %write(handle,force)
end