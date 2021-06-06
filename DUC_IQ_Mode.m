% EXAMPLE FOR IQ MODE HALF IN PROTEUS USING VISA
%===================================================
% This example sets the IQ mode for the designated channel and downloads a
% complex(IQ) waveform to be applied to the built-in IQ modulator in the 
% 'HALF' Mode. This mode offers twice the modulation BW in half of the
% channels.

clc;

fprintf(1, 'INITIALIZING SETTINGS\n');

% Communication Parameters
connStr = '192.168.1.48'; % your IP here
paranoia_level = 0; % 0, 1 or 2

%% Create Administrator
inst = TEProteusInst(connStr, paranoia_level);
fprintf('\n');

res = inst.Connect();
assert (res == true);

% Identify instrument using the standard IEEE-488.2 Command
idnstr = inst.identifyModel();
fprintf('\nConnected to: %s\n', idnstr);

% Get options using the standard IEEE-488.2 Command
optstr = inst.getOptions();

% Get maximum sample rate for target instrument
samplingRate = inst.getMaxSamplingRate();
% Get granularity
granul = inst.getGranularity(idnstr, optstr);

% Carrier Frequency
cfr = 1.0E+09; 
% If sampling rate lower than 2.5GHz, NCO frequency set to 500MHz
if samplingRate < 2.5E+9
    ncoFreq = 500E+6;
end
% Set offset to any positive or negative frequency to shift carrier
fOffset = 0.0;
% Set initial phase for NCO
phase = 0.0 ;
%Set Target Channel
channel = 1;
%Set Target Segment
segment = 1;
% select reversed spectrum for generation in second Nyquist Zone
reverse = false;
% Boost Output Power by 6dB
apply6db = true;

%Wfm Calculation
fprintf(1, 'Calculating WAVEFORM\n');

%ANALOG & DIGITAL MODULATION SETTINGS
% modType Analog:
% 0             AM
% 1             FM
% 2             PM
% 3             SSB; 
% 4             CHIRP;
% modType Digital: 
% 5             QPSK
% 6             QAM16
% 7             QAM32
% 8             QAM64
% 9             QAM128
%10             QAM256
%11             QAM512
%12             QAM1024
modType = 10; 

% See CalculateAnalogModWfm Function to know the meaning of the param1 and
% param2 variables depending on the modulation scheme.
param1 = 1E6; %Peak Frequency Deviaton in Hz
param2 = 1E5; %Modulaiton frequency in HZ
minCycles = 11; %Prime number is better
% Parameters for Digital Modulation
numOfSymbols = 1024;
symbolRate = 30E6;
rollOff = 0.35;
% Interpolation according to DUC interpolation factor
interpol = 4; %8X interpolation factor
% Waveform granularity applies to the combined I/Q waveform so actual
% granularity for each component is granul / 2.
gCorr = 1;
intCorr = 1;

if modType <= 3
    % Calculate analog modulation I/Q waveforms
    wfmIq = CalculateAnalogModWfm(  modType,...
                                    minCycles,...
                                    samplingRate,...
                                    interpol / intCorr,...
                                    granul / gCorr,...
                                    param1,...
                                    param2);
else
    % Calculate QPSK/QAM I/Q Waveforms
    wfmIq = CalculateDigitalModWfm( modType,...
                                    numOfSymbols,...                        
                                    symbolRate,...
                                    rollOff,... 
                                    samplingRate,... 
                                    interpol / intCorr);
    % wfmIq length is not adjusted for granularity to optimize accuracy for
    % symbol rate. Howeverm symbol rate must be adjusted for signal loop
    % consistency. Actual Symbol Rate must be calculated and used in 
    % analysis of the signal.                          
    actualSymbR = samplingRate / interpol * numOfSymbols / length(wfmIq);
    fprintf('\nActual Symbol Rate = to: %d\n', actualSymbR);
end

% I and Q waveforms
myWfmI = real(wfmIq);
myWfmQ = imag(wfmIq);
clear wfmIq;
% Negative Q waveform for inverse spectrum
if reverse
    myWfmQ = -myWfmQ;
end

% Frequency Offset applied
[myWfmI,  myWfmQ] = ApplyFreqOffset(    fOffset,...
                                        samplingRate / interpol,...
                                        myWfmI,...
                                        myWfmQ);

% I/Q data interleaving to a single array for downloadg
fprintf(1, 'I/Q INTERLEAVING\n');
% Envelope normalization to avoid DAC clipping
[myWfmI,  myWfmQ] = NormalIq2(myWfmI, myWfmQ);  

% If necessary, wfm repetated for waveform granularity
myWfmI = trimGran(myWfmI, granul);
myWfmQ = trimGran(myWfmQ, granul);

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Reset AWG
inst.SendCmd('*CLS');
inst.SendCmd('*RST');

% Set sampling rate for AWG to maximum.
inst.SendCmd([':FREQ:RAST ' num2str(2.5E9)]);

% The Half mode requires setting two channels
inst.SendCmd(sprintf(':INST:CHAN %d', channel));
% Interpolation factor for I/Q waveforms set to X4
inst.SendCmd(':SOUR:INT X4');
inst.SendCmd(':MODE DUC');
inst.SendCmd(':IQM HALF');

inst.SendCmd(sprintf(':INST:CHAN %d', channel + 1));
% Interpolation factor for I/Q waveforms set to X4
inst.SendCmd(':SOUR:INT X4');
inst.SendCmd(':MODE DUC');
inst.SendCmd(':IQM HALF');

inst.SendCmd([':FREQ:RAST ' num2str(samplingRate)]);
% DAC Mode set to 'DUC' and IQ Modulation mode set to 'ONE'
 

% Waveform Downloading
% *******************
inst.SendCmd(':TRAC:DEL:ALL');
fprintf(1, 'DOWNLOADING WAVEFORM\n');
res = SendWfmToProteus(inst, channel, segment, myWfmI, 16);
res = SendWfmToProteus(inst, channel + 1, segment + 1, myWfmQ, 16);
fprintf(1, 'WAVEFORM DOWNLOADED!\n');
clear myWfmI myWfmQ;

% Select segment for generation
fprintf(1, 'SETTING AWG OUTPUT\n');
% Q Channel
inst.SendCmd(sprintf(':INST:CHAN %d', channel + 1));
inst.SendCmd(sprintf(':FUNC:MODE:SEGM %d', segment + 1))
% NCO frequency and phase setting
inst.SendCmd(sprintf(':NCO:CFR2 %d', cfr));
inst.SendCmd(sprintf(':NCO:PHAS2 %d', phase));
if apply6db
    inst.SendCmd(':NCO:SIXD2 ON');   
else
    inst.SendCmd(':NCO:SIXD2 OFF');    
end

% I Channel
inst.SendCmd(sprintf(':INST:CHAN %d', channel));
inst.SendCmd(sprintf(':FUNC:MODE:SEGM %d', segment))
% NCO frequency and phase setting
inst.SendCmd(sprintf(':NCO:CFR2 %d', cfr));
inst.SendCmd(sprintf(':NCO:PHAS2 %d', phase));
if apply6db
    inst.SendCmd(':NCO:SIXD2 ON');   
else
    inst.SendCmd(':NCO:SIXD2 OFF');    
end

% Output volatge set to MAX
inst.SendCmd(':SOUR:VOLT MAX');
% Activate outpurt and start generation
inst.SendCmd(':OUTP ON');

fprintf(1, 'SETTING SAMPLING CLOCK\n');
% Set sampling rate for AWG as defined in the preamble.
inst.SendCmd([':FREQ:RAST ' num2str(samplingRate)]);


% It is recommended to disconnect from instrumet at the end
inst.Disconnect();    
clear inst;
clear;
fprintf(1, 'END\n');

function finalWfm = trimGran(inWfm, granularity)
    % trimGran - Adjust wfm length for granularity
    %
    % Synopsis
    %   finalWfm = trimGran(inWfm, granularity)
    %
    % Description
    %   Repeat waveform the minmum number of times to meet the
    %   waveform length granularity criteria
    %
    % Inputs ([]s are optional)
    %   (double) inWfm  Input waveform
    %   (int16)  granularity
    %
    % Outputs ([]s are optional)
    %   (double) finalWfm Adjusted waveform

    baseL = length(inWfm);
    finaL = lcm(baseL, granularity);
    
    finalWfm = zeros(1, finaL);
    pointer = 1;
    
    while pointer < finaL
        finalWfm(pointer : (pointer + baseL -1)) = inWfm;
        pointer = pointer + baseL;        
    end

end

function [rotI,  rotQ] = ApplyFreqOffset(fOffset, sampleRate, wfmI, wfmQ)    
    
    wfmL = length(wfmI);
    fRes = sampleRate / wfmL;
    fOffset = round(fOffset / fRes) * fRes;

    cplexWfm = wfmI + 1i * wfmQ;
    clear wfmI wfmQ;
    angleArray = 0:(wfmL - 1);
    angleArray = 2 * pi * fOffset * angleArray;
    angleArray = angleArray / sampleRate;
    
    angleArray = exp(1i * angleArray);
    
    cplexWfm = cplexWfm .* angleArray;
    clear angleArray;
    
    rotI = real(cplexWfm);
    rotQ = imag(cplexWfm);    
end

function [normI,  normQ] = NormalIq2(wfmI, wfmQ)    
    
    maxPwr = max(abs(wfmI));
    
    if maxPwr < max(abs(wfmQ))
        maxPwr = max(abs(wfmQ));
    end
    
    normI = wfmI / maxPwr;
    normQ = wfmQ / maxPwr;
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
    myWfm = instHandle.Quantization(myWfm, dacRes);
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';
    
    if dacRes == 16
        instHandle.SendBinaryData(prefix, myWfm, 'uint16');
    else
        instHandle.SendBinaryData(prefix, myWfm, 'uint8');
    end   
    
    result = length(myWfm);
end

function waveform = CalculateAnalogModWfm(  modType,...
                                            minCycles,...
                                            sampleRate,...
                                            interpol,...
                                            granul,...
                                            param1,...
                                            param2)
   
    %ANALOG MODULATION WAVEFORM CALCULATION
    % modType = 0, AM; 1, FM; 2, PM; 3, SSB;    

    %AM SETTINGS
    amModIndex = param1; %Modulation Index in %
    amModFreq = param2; %Modulation frequency in HZ

    %FM SETTINGS
    fmFreqDev = param1; %Peak Frequency Deviaton in Hz
    fmModFreq = param2; %Modulaition frequency in HZ

    %PM SETTINGS
    pmPhaseDev = param1; %Peak Phase Deviaton in Rads
    pmModFreq = param2; %Modulation frequency in HZ

    %SSB SETTINGS
    ssbModFreq = param2; %Modulation frequency in HZ

    %CHIRP SETTINGS
    chirpSweepRange = param1;
    chirpSweepTime = param2;

    %Waveform Length Calculation    
    modFreq = amModFreq;

    if modType == 1
        modFreq = fmModFreq;
    elseif modType == 2
        modFreq = pmModFreq;
    elseif modType == 3
        modFreq = ssbModFreq;
    elseif modType == 4
        modFreq = 1 / chirpSweepTime;
    end

    actualSR = sampleRate / interpol;
    if modType ~= 4 
        numOfSamples = round(actualSR / abs(modFreq / minCycles));
    else
        numOfSamples = round(actualSR / abs(modFreq));
    end
    totalNumOfSamples = numOfSamples;
    
    % As samples sent to the instrument are twice the number of complex
    % samples, granul must be defined as half the actual number
    
    while modType ~= 4 && mod(totalNumOfSamples, granul) ~= 0
        totalNumOfSamples = totalNumOfSamples + numOfSamples;
    end

    numOfSamples = totalNumOfSamples;
    fRes = actualSR / numOfSamples;

    % Round modFreq to the nearest integer number of Cycles

    modFreq = round(modFreq / fRes) * fRes;

    %Waveform calculation
    fprintf(1, 'WAVEFORM CALCULATION\n');

    waveform = 0: (numOfSamples - 1);
    waveform = (1 / actualSR) .* waveform;
    waveform = waveform - (numOfSamples / (2 * actualSR));

    if modType == 0
        waveform = 1 + amModIndex/100 .* sin(2 * pi * modFreq * waveform);
        waveform = waveform + 1i * waveform;
    elseif modType == 1
        fmFreqDev = round(fmFreqDev / fRes) * fRes;
        freqInst = fmFreqDev / modFreq * sin(2 * pi * modFreq * waveform);
        waveform = cos(freqInst) + 1i * sin(freqInst); 
        clear freqInst;
    elseif modType == 2
        phaseInst = pmPhaseDev * sin(2 * pi * modFreq * waveform);
        waveform = cos(phaseInst) + 1i * sin(phaseInst);
        clear phaseInst;
    elseif modType == 3
        waveform = 2 * pi * modFreq * waveform;
        waveform = cos(waveform) + 1i * sin(waveform);
    elseif modType == 4
        chirpSweepRange = chirpSweepRange / 2;
        chirpSweepRange = round(chirpSweepRange / fRes) * fRes;
        freqInst = (actualSR * chirpSweepRange / numOfSamples) * waveform;
        freqInst = 2 * pi * freqInst .* waveform;
        waveform = cos(freqInst) + 1i * sin(freqInst);
        clear freqInst;    
        waveform = trimGran(waveform, granul);    
    end

    % waveform conditioning:    
    waveform = waveform./((mean(abs(waveform).^2))^0.5);

end
function [dataOut] = CalculateDigitalModWfm(    modType,...
                                                numOfSymbols,...
                                                symbolRate,...
                                                rollOff,... 
                                                sampleRate,... 
                                                interpol)
                            
    % modType     Modulation
    % 5             QPSK
    % 6             QAM16
    % 7             QAM32
    % 8             QAM64
    % 9             QAM128
    %10             QAM256
    %11             QAM512
    %12             QAM1024
    
    if modType == 5
        bitsPerSymbol = 2;
    elseif modType == 6
        bitsPerSymbol = 4;
    elseif modType == 7
        bitsPerSymbol = 5;
    elseif modType == 8
        bitsPerSymbol = 6;
    elseif modType == 9
        bitsPerSymbol = 7;
    elseif modType == 10
        bitsPerSymbol = 8;
    elseif modType == 11
        bitsPerSymbol = 9;
    elseif modType == 12
        bitsPerSymbol = 10;
    else
        bitsPerSymbol = 2;
    end
                                
    % Waveform Length Calculation
    sampleRate = sampleRate / interpol;
       
    % Create IQ for QPSK/QAM    
    % accuracy is the length of the shaping filter
    accuracy = 32;
    fType = 'sqrt'; % 'normal' or 'sqrt'
    % Get symbols in the range 1..2^bps-1
    data = getRnData(numOfSymbols, bitsPerSymbol);
    % Map symbols to I/Q constellation locations
    [dataI, dataQ] = getIqMap(data, bitsPerSymbol);
    % Adapt I/Q sample rate to the AWG's
    oversampling = round(sampleRate / symbolRate);
    dataI = expanData(dataI, oversampling);
    dataQ = expanData(dataQ, oversampling);
    % Calculate baseband shaping filter
    rsFilter = rcosdesign(rollOff,accuracy,oversampling, fType);
    % Apply filter through circular convolution
    dataI = cconv(dataI, rsFilter, length(dataI));
    dataQ = cconv(dataQ, rsFilter, length(dataQ));
    % Output waveforfm must be made of complex samples
    dataOut = dataI + 1i * dataQ;
end

function dataOut = getRnData(nOfS, bPerS)

    maxVal = 2 ^ bPerS;
    dataOut = maxVal * rand(1, nOfS);
    dataOut = floor(dataOut);
    dataOut(dataOut >= maxVal) = maxVal - 1;    
end

function [symbI, symbQ] = getIqMap(data, bPerS)
   
    if bPerS == 5 % QAM32 mapping
        lev = 6;
        data = data + 1;
        data(data > 4) = data(data > 4) + 1;
        data(data > 29) = data(data > 29) + 1; 
    
    elseif bPerS == 7 % QAM128 mapping      
        lev = 12;
        data = data + 2;
        data(data > 9) = data(data > 9) + 4;
        data(data > 21) = data(data > 21) + 2;
        data(data > 119) = data(data > 119) + 2;
        data(data > 129) = data(data > 129) + 4;
        
     elseif bPerS == 9 % QAM512 mapping       
        lev = 24;
        data = data + 4;
        data(data > 19) = data(data > 19) + 8;
        data(data > 43) = data(data > 43) + 8;
        data(data > 67) = data(data > 67) + 8;
        data(data > 91) = data(data > 91) + 4;
        data(data > 479) = data(data > 479) + 4;
        data(data > 499) = data(data > 499) + 8;
        data(data > 523) = data(data > 523) + 8;
        data(data > 547) = data(data > 547) + 8;            
    else
        lev = 2 ^ (bPerS / 2); % QPSK, QAM16, QAM64, QAM256, QAM1024      
    end
    
    symbI = floor(data / lev);
    symbQ = mod(data, lev);
    lev = lev / 2 - 0.5;   
    symbI = (symbI - lev) / lev;
    symbQ = (symbQ - lev) / lev;


end

function dataOut = expanData(inputWfm, oversampling)

    dataOut = zeros(1, oversampling * length(inputWfm));
    dataOut(1:oversampling:length(dataOut)) = inputWfm;

end

    