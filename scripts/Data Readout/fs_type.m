function o = fs_type(type)
% Returns the type as parsed in FS file types, as a string that can be used
% for the matlab function "typecast"
%
% Usage:
% o = fs_type(type);

% File types
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;

if(type < 1 || type > 8 || type == FS_CHAR)
	o = 'char';
elseif(type == FS_UCHAR)
	o = 'int8';
elseif(type == FS_INT)
	o = 'int32';
elseif(type == FS_UINT)
	o = 'uint32';
elseif(type == FS_FLOAT)
	o = 'float';
elseif(type == FS_DOUBLE)
	o = 'double';
elseif(type == FS_INT64)
	o = 'int64';
elseif(type == FS_UINT64)
	o = 'uint64';
end