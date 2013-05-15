%get a device
matlabrc;
h = haptikdevice;

%run simulation for 10 seconds
tic
while toc < 10
    
    %read probe position
    pos = read_position(h);


    %check for collision and send back force feedback
    if pos(2)<0
        write(h, -1 * [0 pos(2) 0]);
    else
        write(h,[0 0 0]);
    end
    
end

close(h);
clear h;
clear pos;

