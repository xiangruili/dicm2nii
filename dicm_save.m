function dicm_save(img, fname, s)
% msg = DICM_SAVE(img, dicomFileName, info_struct);
% Save img into dicom file, using tags stored in struct. The dicom image type is
% based on the data type of img.
%
% The current version overrides TransferSyntaxUID with '1.2.840.10008.1.2.1':
% Little Endian, Explicit VR, no compression, and for tags in dicm_dict.m only. 
%
% See also DICM_DICT, DICM_HDR, DICM_IMG

% 200109 (yymmdd) Write it ( xiangrui.li at gmail.com)

persistent dict;
narginchk(2, 3);
if ~isnumeric(img), error('Provide img array as 1st input.'); end
if ~ischar(fname) && ~isstring(fname), error('Need a string as file name.'); end

fid = fopen(fname, 'w', 'l');
if fid<0, error('Failed to open file %s', fname); end
closeFile = onCleanup(@() fclose(fid)); % auto close when done or error
if isempty(dict)
    try vndor = s.Manufacturer; catch, vndor = ''; end
    dict = dicm_dict(vndor); % no FileMetaInformationGroupLength
end

s.FileMetaInformationVersion = [0 1]';
s.TransferSyntaxUID = '1.2.840.10008.1.2.1'; % add or overwrite
if ~isfield(s, 'SOPInstanceUID'), s.SOPInstanceUID = dicm_uid(s); end
s.MediaStorageSOPInstanceUID = s.SOPInstanceUID;

sz = size(img); sz(numel(sz)+1:4) = 1;
if sz(4)>1, s.NumberOfFrames = sz(4); end
s.Rows = sz(1);
s.Columns = sz(2);
typ = class(img);
switch typ % only first 2 supported by dicom standard
    case {'int8'  'uint8' }, vr = 'OB'; nBit = 8; 
    case {'int16' 'uint16'}, vr = 'OW'; nBit = 16;
    case {'int32' 'uint32'}, vr = 'OL'; nBit = 32;
    case {'int64' 'uint64'}, vr = 'OV'; nBit = 64;
    case 'single',           vr = 'OF'; nBit = 32;
    case 'double',           vr = 'OD'; nBit = 64;
    otherwise, error('Unknown image type: %s', typ);
end
if strfind(typ, 'int'), s.PixelRepresentation = typ(1)~='u'; end %#ok
s.BitsAllocated = nBit;
s.BitsStored = nBit;
s.HighBit = nBit - 1;
s.LargestImagePixelValue = max(img(:));
s.SmallestImagePixelValue = min(img(:));

fwrite(fid, zeros(1,16), 'int64'); % 128 byte
fwrite(fid, 'DICM', 'char*1'); % signature
fwrite(fid, [0 0 0], 'uint32'); % 12-byte for g2Len at 132
g2End = write_meta(fid, s, dict);
write_tag(fid, [32736 16], vr, permute(img, [2 1 3:ndims(img)])); % PixelData
fseek(fid, 132, 'bof'); % very first tag
write_tag(fid, [2 0], 'UL', g2End-144); % FileMetaInformationGroupLength

%% Write most meta info in struct s
function g2End = write_meta(fid, s, dict)
g2End = 0; % only used by main func for once
flds = fieldnames(s);
flds(strcmp(flds, 'PixelData')) = [];
N = numel(flds);
ind = zeros(1, N);
for i = 1:N
    j = find(strcmp(flds{i}, dict.name));
    if isempty(j), continue; end
    if numel(j) > 1
        pub = mod(dict.group(j), 2) == 0;
        if any(pub), j = j(pub); end
    end
    ind(i) = j(1);
end
for i = sort(ind(ind>0)) % in the order of tags
    if g2End<1 && dict.group(i)>2, g2End = ftell(fid); end
    write_tag(fid, [dict.group(i) dict.element(i)], dict.vr{i}, s.(dict.name{i}), dict);
end

%% Write a tag
function write_tag(fid, tag, vr, val, dict)
fwrite(fid, tag, 'uint16');
fwrite(fid, vr, 'char*1');
if strcmp(vr, 'SQ')
    fwrite(fid, [0 65535 65535], 'uint16'); % FFFF FFFF, lazy SQ length
    if isstruct(val)
        if isfield(val, 'Item_1'), val = struct2cell(val);
        else, val = {val};
        end
    end
    for i = 1:numel(val)
        if ~isstruct(val{i}), continue; end
        fwrite(fid, [65534 57344 65535 65535], 'uint16'); % FFFE E000 (Item) & len
        write_meta(fid, val{i}, dict); % recuisive call
        fwrite(fid, [65534 57357 0 0], 'uint16'); % FFFE E00D, ItemDelimitationItem
    end
    fwrite(fid, [65534 57565 0 0], 'uint16'); % FFFE E0DD, SequenceDelimitationItem
    return;
end

if iscellstr(val), val = sprintf('%s\n', val{:}); end %#ok
if strcmp(vr, 'DS') || strcmp(vr, 'IS')
    fmt = 'char*1';
    val = sprintf('%.16g\\', val); val(end) = '';
    n = numel(val);
else
    [fmt, bpv] = vr2fmt(vr);
    n = numel(val) * bpv;
end

len16 = 'AE AS AT CS DA DS DT FD FL IS LO LT PN SH SL SS ST TM UI UL US';
if isempty(strfind(len16, vr)) %#ok<*STREMP>
    fmtLen = 'uint32';
    fwrite(fid, 0, 'uint16');
else
    fmtLen = 'uint16';
end
odd = mod(n, 2);
fwrite(fid, n+odd, fmtLen); % length is even
fwrite(fid, val, fmt);
if odd, fwrite(fid, 0, 'uint8'); end

%% Return format str and Bytes per value from VR
function [fmt, bpv] = vr2fmt(vr)
switch vr
    case {'AE' 'AS' 'CS' 'DA' 'DT' 'LO' 'LT' 'PN' 'SH' 'ST' 'TM' 'UI' 'UT'}
               bpv = 1; fmt = 'char*1';
    case 'US', bpv = 2; fmt = 'uint16';
    case 'OB', bpv = 1; fmt = 'uint8';
    case 'FD', bpv = 8; fmt = 'double';
    case 'SS', bpv = 2; fmt = 'int16';
    case 'UL', bpv = 4; fmt = 'uint32';
    case 'SL', bpv = 4; fmt = 'int32';
    case 'FL', bpv = 4; fmt = 'single';
    case 'AT', bpv = 2; fmt = 'uint16';
    case 'OW', bpv = 2; fmt = 'uint16';
    case 'OF', bpv = 4; fmt = 'single';
    case 'OD', bpv = 8; fmt = 'double';
    otherwise, bpv = 1; fmt = 'uint8';
end

%% Generate 'unique' dicom ID
function uid = dicm_uid(s)
try pre = s.SOPInstanceUID(1:27);
catch, pre = '1.3.12.2.1107.5.2.43.167002'; % CCBBI Prisma
end
uid = sprintf('%s.1.%s%07.0f', pre, datestr(now, 'yyyymmddHHMMSSfff'), rand*1e7);
