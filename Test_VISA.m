% BASIC EXAMPLE FOR CONNECTION TO PROTEUS USING VISA
%===================================================
% VISA Communications from MATLAB requires the Instrument Control Toolbox

clc;

% Define IP Address for Target Proteus device descriptor
% VISA "Socket-Based" TCP-IP Device. Socket# = 5025
connStr = '192.168.1.48';
paranoia_level = 0;

% Connect to Instrument and get visa Handle
% Using the TEProteusInst library
inst = TEProteusInst(connStr, paranoia_level);
res = inst.Connect();
assert (res == true);

% Identify instrument using the standard IEEE-488.2 Command
idnstr = inst.identifyModel();
fprintf('\nConnected to: %s\n', idnstr);

% Get options using the standard IEEE-488.2 Command
optstr = inst.getOptions();

for i=1:length(optstr)
    fprintf('\nOption #%d Installed: %s', i, char(optstr(i)));     
end

% Get Number of Channels
numchan = inst.getNumOfChannels(idnstr);
fprintf('\n\nNumber of Channels = %d\n', numchan);

% Get min sampling rate
minsr = inst.getMinSamplingRate();
fprintf('\nMinimum Sample Rate = %d samples/sec\n', minsr);

% Get max sampling rate
maxsr = inst.getMaxSamplingRate();
fprintf('\nMaximum Sample Rate = %d samples/sec\n', maxsr);

% Get granularity
granul = inst.getGranularity(idnstr, optstr);
fprintf('\nGranularity = %d samples\n', granul);

% Get Sample Resolution
dacres = inst.getDacResolution();
fprintf('\nSample Resolution = %d bits\n', dacres);

% Disconnect, close VISA handle and destroy handle  
inst.delete();
