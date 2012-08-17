function v12 = interp_vec(v1, v2)
v1l = length(v1);
v2l = length(v2);

if(v1l > v2l)
	v1(v2l:end) = [];
elseif(v2l > v1l)
	v2(v1l:end) = [];
end

vl = length(v1);

v12 = zeros(vl*2, 1);

v12(1:2:end) = v1;
v12(2:2:end) = v2;