
function [dr,frsize,usercv,rhnframes,nslices,hdr] = rawloadHD_si(fn,frsize,ep,ishosp)


% 1 	version
% 40 	date, time, etc
% 88 	hdr
% 19  	usercv
% 439 	hdr 
% 11	garbage
% 58324 garbage
% 1044 	garbage


if (nargin < 4) | (length(ishosp)==0)
	ishosp=0;
end;

%
%  **** For more information, the P-file header is well documented in 
%	the LX ESE users manual.
%

if (nargin < 1)	% Get last p-file if function exists.
	if (exist('lastpfile'))
	    fn = lastpfile;
	    tt = sprintf('No file given... using most recent file, %s.',fn);
	    %disp(tt);
	else
		error('You must specify a P-file');
	end;
end;

% default to normal precision
if (nargin < 3) | (length(ep)==0)
  ep = 0;
end

fip = fopen(fn,'r'); % Took away 'b' with Excite upgrade...
if fip == -1
  tt = sprintf('File %s not found\n',fn);
  %dispsafe(tt);
  return;
end

ver = fread(fip,1,'float');
%dispsafe('Raw file version:');
 
if(ver==11) 			tothdrsize = 66072;  % 60464 in EX  %39940 in LX;
elseif(14<=ver && ver <=16)	tothdrsize = 145908; % 14.x
elseif(20<=ver && ver <21)	tothdrsize = 149788; % 20.x
end
 

hdrrecsize = 2048;
daqtabsize = 20480;
if ((ver > 7) & (ishosp==0))
	daqtabsize = 2*20480;	% Fixed?!
end;
othhdrsize = 4096+4096+daqtabsize+2052+2052+2048+1000;  % 1000 is empirical
examhdrsize = 1036;	% Was 1040, changed Oct 27/03 
serieshdrsize = 984; %1028;  Could this be the 44-byte difference??
imagehdrsize = 1044;

% These should add to tothdrsize.
hdrrecsize+othhdrsize+examhdrsize+serieshdrsize+imagehdrsize;

% ===========================================================
% Extract Date and Time from header.
% ===========================================================

nstr = 40;
%nstr = 44;

hstr = fread(fip,nstr,'char');
hstr = char(hstr(:))';
%dispsafe(hstr);

dstr = sprintf('%c',hstr(13:22));
if (dstr(7:8)=='10')
	dstr = [dstr(1:6) '200' dstr(9)];
end;
tstr = sprintf('%c',hstr(23:30));

%disp(' ');
%disp('============================================================');
%disp('	RawloadEX --> EXCITE Version.				');
%disp('============================================================');
dt = sprintf('Raw File:  %s ',fn);
%dispsafe(dt);
dt = sprintf('Data Acquisition Date/Time: %s, %s',dstr,tstr);
%dispsafe(dt);


% ===========================================================
% Read part of the header up to and including rhuser0
% ===========================================================
hdr = fread(fip,88,'int16');

% for n = 1: 88
%    disp(sprintf('hdr(%d) = %d', n, hdr(n)));
% end


if (nargin < 2) | (length(frsize)==0)
  frsize = hdr(19);
  dt = sprintf('Detected Framesize of %d',frsize);
else
  dt = sprintf('Specified Framesize of %d',frsize);
end
%dispsafe(dt);
if (nargin < 4) | (length(ishosp)==0)
	ishosp=0;
end;

rhnframes=hdr(16);
nslices=hdr(13);
nechoes=hdr(14);
nex=hdr(15);

dt = sprintf('%d frames, %d slice(s), %d echo(es), %d NEX',rhnframes,nslices,nechoes,nex);
%dispsafe(dt);

rawsize = 32768-sign(hdr(37))*32768 + hdr(37)+65536*hdr(38);
	% rawsize is a 32-bit value, but stored as 2 16-bit words.
	% The above is the conversion.
bytes_per_sample = 4;
nframes = floor(rawsize / frsize / bytes_per_sample +0.5);

dt = sprintf('%d Total Frames expected',nframes);
% dispsafe(dt);
 
clear hdr;


% ===========================================================
% Read rhuser1-rhuser19
% ===========================================================
usercvsize = 19*4;
usercv = fread(fip,usercvsize/4,'float');

%for n = 1: 19
%   disp(sprintf('usercv(%d) = %f', n, usercv(n)));
%end

% ===========================================================
% Read rest of RDB header.
% ===========================================================
hdr = fread(fip,(hdrrecsize-nstr-usercvsize-176)/4,'int32'); % (hdrrecsize-nstr-usercvsize-176)/4 = 439

%for n = 1: (hdrrecsize-nstr-usercvsize-176)/4
%   disp(sprintf('hdr(%d) = %d', n, hdr(n)));
%end


prescanpars = hdr(30:37);	
tt = sprintf('Prescan Parameters:   AX=%ld  TG=%d  R1=%d  R2=%d',prescanpars(8),prescanpars(7),prescanpars(5),prescanpars(6));
%dispsafe(tt);

fread(fip,11,'float');

% ===========================================================
% Read Next Few headers.
% ===========================================================
hdr = fread(fip,othhdrsize+examhdrsize+serieshdrsize,'char'); % othhdrsize+examhdrsize+serieshdrsize = 58324

% ===========================================================
% Read Image Header, and interpret/display some Info.
% ===========================================================

hdr = fread(fip,floor((imagehdrsize)),'uint8'); % imagehdrsize = 1044

% ======= TR, TE =========
tr = bytes2int32(hdr(201:204));
te = bytes2int32(hdr(209:212));
tt = sprintf('TR=%8.3f ms,  TE=%8.3f ms',tr/1000,te/1000);
%dispsafe(tt);	 TR/TE NOT WORKING IN EXCITE HEADER YET!!

% ====== PSD Name ========
tt = sprintf('PSD = %s',char(hdr(321:353))');
%dispsafe(tt);	PSD NOT WORKING WITH EXCITE HEADER YET!!


% ====== Trigger Delay =======
tdel = bytes2int32(hdr(233:236));
hr = bytes2int16(hdr(231:233));
tt = sprintf('Cardiac Info:  Rate = %8.3f bpm  Delay = %8.3f ms',hr,tdel/1000);
%dispsafe(tt);	CARDIAC INFO NOT WORKING WITH EXCITE HEADER YET!!
%dispsafe('============================================================');


% ===========================================================
% Read all the data
% ===========================================================


fseek(fip,tothdrsize,-1);
%if ep == 1 
%  dr = fread(fip,nframes*frsize*2,'int32');
%  dd = fread(fip,inf,'int32');
%else
%  dr = fread(fip,nframes*frsize*2,'int16');
%  dd = fread(fip,inf,'int16');
%end
if ep == 1 
  dr = fread(fip,inf,'int32');
else
  dr = fread(fip,inf,'int16');
end


ll = length(dr);
dr = dr(1:2:ll)+i*dr(2:2:ll);
nframes = floor(ll/2/frsize); 
%ld = length(dd);
%tt = sprintf('%d complex samples discarded at end of file',ld/2);
%dispsafe(tt);

dr = dr(1:nframes*frsize); % redundant
dr = reshape(dr,frsize,nframes);

fclose(fip);

return;





function intout = bytes2int32(bytes)
%
%	Convert 4 bytes to an int32
%
intout = (((bytes(1)*256)+bytes(2)*256)+bytes(3)*256)+bytes(4);
if (intout >= 32768*65536)
	intout = intout-65536*65536;
end;



function intout = bytes2int16(bytes)
%
%	Convert 4 bytes to an int32
%
intout = (bytes(1)*256)+bytes(2);
if (intout >= 32768)
	intout = intout-65536;
end;



