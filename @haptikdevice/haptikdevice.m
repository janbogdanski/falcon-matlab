function h = haptikdevice(arg1,arg2,arg3)

% default parameters
 id = -1;
 cb = 0;
 rate = 1000;
 
% input arguments
if nargin >= 1
    if strcmp(class(arg1),'function_handle')
 		cb = arg1;
    else
 		id = arg1;
    end
end
 if nargin >= 2
 	if strcmp(class(arg2),'function_handle')
 		cb = arg2;
    else
 		rate = arg2;
    end
end
if nargin >= 3
	rate = arg3;
end


% create class
h.id = haptik_matlab(1,id);
h.t = 0;
h = class(h,'haptikdevice');

% start device
haptik_matlab(2,h.id);

%if callback then set up timer
if strcmp(class(cb),'function_handle')
	h.t = timer('ExecutionMode','fixedRate','Period', rate/1000000,'TimerFcn',{@internal_haptik_callback,cb,h} );
	start(h.t);
end


% redirector to callback
function internal_haptik_callback(obj,event,cb,h)
feval(cb,h);
