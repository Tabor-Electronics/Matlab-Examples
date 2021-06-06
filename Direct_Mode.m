% EXAMPLE FOR DIRECT MODE
%===================================================
% This example calculates up to 4 different signals and download them into
% each available channel in the target Proteus device.
% 
% The basic waveform is an square waveform using the full DAC range and it
% is downloaded to channel #1. For each channel, the waveform is calculated
% by integration of the previous waveform in a similar way to some analog
% signal generators, where the triangular wave is obtained by integration
% of an square wave, and the sinusoidal waveform is obtained by integration
% of the triangular wave. Channel #4, when available, will generate a
% "cosine" wave obatined by integration of the sinewave assigned to channel
% #3.

clc;

fprintf(1, 'INITIALIZING SETTINGS\n');

% Communication Parameters
connStr = '192.168.1.48'; % your IP here
paranoia_level = 2; % 0, 1 or 2

%% Create Administrator
inst = TEProteusInst(connStr, paranoia_level);
fprintf('\n');

res = inst.Connect();
assert (res == true);

% Identify instrument using the standard IEEE-488.2 Command
idnstr = inst.identifyModel();
fprintf('\nConnected to: %s\n', idnstr);

% Reset AWG
inst.SendCmd('*CLS');
inst.SendCmd('*RST');

% Get options using the standard IEEE-488.2 Command
optstr = inst.getOptions();

% Get granularity
granul = inst.getGranularity(idnstr, optstr);

% Get Number of Channels
numOfChannels = inst.getNumOfChannels(idnstr);
% This example is written to handle a maximum of 4 channels
if numOfChannels > 4
    numOfChannels = 4;
end

% Get maximum sample rate for target instrument
samplingRate = inst.getMaxSamplingRate();
% Sample rate is limited to 2.5GS/s so it is alway possible to generate
% waveforms for all the channels no matter the target Proteus model.
% For the P908X sampling rate can reach 9GS/s in direct mode for all the
% channels. For the P948X models, half the channels are available for
% direct generation over 2.5GS/s without interpolation,
if samplingRate > 2.5E9
    samplingRate = 2.5E9;
end

fprintf(1, 'Calculating WAVEFORMS\n');

minCycles = 1;
period = 1.25E-6;
segment = 1;

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Set sampling rate for AWG to maximum.
inst.SendCmd([':FREQ:RAST ' num2str(samplingRate)]);
%inst.SendCmd(':TRAC:FORM U8');

% Get the default DAC resolution for the current settings.
dacRes = inst.getDacResolution();
% Calculate basic square wave
myWfm = getSquareWfm(   samplingRate,... 
                        minCycles,...
                        period,...
                        granul);

for channel = 1:numOfChannels
    %Select Channel
    inst.SendCmd(sprintf(':INST:CHAN %d', channel));
    % DAC Mode set to 'DIRECT" (Default)
    inst.SendCmd(':MODE DIR');
    
    % Segment # processing
    actualSegment = segment;
    % All Proteus models except the P908X share the same waveform memory
    % bank among channel N+1 and N+2, N=0..NumOfChannels/2. This means that
    % the same segment number cannot be used for this pair of channels. In
    % this case the designated segment is used for the odd numbered
    % channels and the next segment is assigned to the even numbered channel
    % of the same pair. All segments can be deleted just once for each pair
    % of channels.
    if mod((channel - 1), 2) == 1
        actualSegment = segment + 1;
    else
        % All segments deleted for current waveform memory bank
        inst.SendCmd(':TRAC:DEL:ALL');
    end    
    
    % Waveform Downloading
    % *******************
    
    fprintf(1, 'DOWNLOADING WAVEFORM FOR CH%d\n', channel);
    res = SendWfmToProteus(inst, channel, actualSegment, myWfm, dacRes);
    fprintf(1, 'WAVEFORM DOWNLOADED!\n');
    % Select segment for generation
    fprintf(1, 'SETTING AWG OUTPUT\n');
    inst.SendCmd(sprintf(':FUNC:MODE:SEGM %d', actualSegment))
    % Output volatge set to MAX
    inst.SendCmd(':SOUR:VOLT MAX');   
    % Activate outpurt and start generation
    inst.SendCmd(':OUTP ON');
    
    % The new waveform is calculated for the next channel
    if channel < numOfChannels
        % Integration
        myWfm = cumsum(myWfm);
        % DC removal
        myWfm = myWfm - mean(myWfm);
        % Normalization to the -1.0/+1.0 range
        myWfm = myWfm / max(abs(myWfm));
    end
end

% It is recommended to disconnect from instrument at the end
inst.Disconnect();    
clear inst;
clear;
fprintf(1, 'END\n');

function sqrWfm = getSquareWfm( samplingRate,... 
                                numCycles,...
                                period,...
                                granularity)
                            
    wfmLength = round(numCycles * period *samplingRate);
    wfmLength = round(wfmLength / granularity) * granularity;
    
    period = wfmLength / numCycles;    
    sqrWfm = 0:(wfmLength - 1);    
    sqrWfm = square(sqrWfm * 2 * pi / period);    
                            
end

function retval = myQuantization (myArray, dacRes)
  
  minLevel = 0;
  maxLevel = 2 ^ dacRes - 1;  
  numOfLevels = maxLevel - minLevel + 1;
  
  retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
  retval = retval + minLevel;
  
  retval(retval > maxLevel) = maxLevel;
  retval(retval < minLevel) = minLevel;

end

function result = SendWfmToProteus( instHandle,...
                                    channel,...
                                    segment,...
                                    myWfm,...
                                    dacRes)

    %Select Channel
    instHandle.SendCmd(sprintf(':INST:CHAN %d', channel));    
    instHandle.SendCmd(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)));        
    % select segmen as the the programmable segment
    instHandle.SendCmd(sprintf(':TRAC:SEL %d', segment));

    % format Wfm
    myWfm = myQuantization(myWfm, dacRes);
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';
    
    if dacRes == 16
        instHandle.SendBinaryData(prefix, myWfm, 'uint16');
    else
        instHandle.SendBinaryData(prefix, myWfm, 'uint8');
    end   
    
    result = length(myWfm);
end  