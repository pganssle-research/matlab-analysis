function hash = basic_hash(string)

hash = 0;
nl = 16*(floor(length(string))+1);
nstring = ones(1, nl);

nstring(1:length(string)) = string;

for i = 1:length(nstring)
    hash = mod((37*(hash + nstring(i))), hex2dec('ffffffffffffffff'));
end


