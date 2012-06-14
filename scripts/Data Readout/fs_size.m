function o = fs_size(type)
% Gives the sizes of various types as per the FS file spec (MCD, PP), etc.
%
% Usage:
% o = fs_size(type);

% File types
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;

if(type < 1 || type > 8)
	o = 1;
elseif(type == FS_CHAR || type == FS_UCHAR)
	o = 1;
elseif(type == FS_INT || type == FS_UINT || type == FS_FLOAT)
	o = 4;
elseif(type == FS_DOUBLE || type == FS_INT64 || type == FS_UINT64)
	o = 8;
end