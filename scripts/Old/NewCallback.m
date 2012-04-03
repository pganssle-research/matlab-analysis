function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
ind = get(event_obj,'DataIndex');

save('eventobj.mat', 'event_obj');

oldval = get(event_obj.Target, 'UserData');
so = size(oldval);
set(event_obj.Target, 'UserData', pos);

output_txt = {['X: ',num2str(pos(1),6)], ...
        ['Y: ',num2str(pos(2),6)],...
        ['Index: ', num2str(ind, 0)]};
    
if(so(1) ~= 0)
	output_txt = {output_txt{:}, ['/\X: ', num2str((oldval(1)-pos(1)), 6)], ['/\Y: ', num2str((oldval(2)-pos(2)), 6)]};
end       

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
