% EXAMPLE FOR IQ MODE TWO IN PROTEUS USING VISA
%==============================================
% This example sets the IQ mode for the designated channel and downloads
% two complex(IQ) waveforms to be applied to the two built-in IQ modulators 
% per channel in the 'TWO' Mode.

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

% Get options using the standard IEEE-488.2 Command
optstr = inst.getOptions();

% Get maximum sample rate for target instrument
samplingRate = inst.getMaxSamplingRate();
% Get granularity
granul = inst.getGranularity(idnstr, optstr);

% Carrier Frequency
ncoFreq1 = 100E6;
ncoFreq2 = 2.0E+09; 
% If sampling rate lower than 2.5GHz, NCO frequency set to 500MHz
%samplingRate = 2.5E+9;
if samplingRate > 5000E+6
    samplingRate = 5000E+6;
end
if samplingRate < 2.5E+9
    ncoFreq2 = 500E+6;
end
% Set offset to any positive or negative frequency to shift carrier
fOffset1 = 0.0E6;
fOffset2 = 0.0E6;

% Set initial phase for NCO
phase1 = 0.0;
phase2 = 0.0;
%Set Target Channel
channel = 1;
%Set Target Segment
segment = 1;
% select reversed spectrum for generation in second Nyquist Zone
reverse = false;
% Boost Output Power by 6dB
apply6db = true;

%Wfm Calculation
fprintf(1, 'Calculating WAVEFORMS\n');

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
modType1 = 8;
modType2 = 5;

minCycles = 11; %Prime number is better
interpol = 8; 
% Effective granularity for waveforms is 1/4 of the waveform memory
% granularity as 4 interleaved samples are stores for each sample time.
gCorrection = 4;
intCorrection = 1;

% See CalculateAnalogModWfm Function to know the meaning of the param1 and
% param2 variables depending on the modulation scheme.
param1_1 = 1E7; %Peak Frequency Deviaton in Hz
param2_1 = 1E5; %Modulation frequency in HZ

numOfSymbols1 = 1024;
symbolRate1 = 15E6;
rollOff1 = 0.35;

% Parameters for second DUC
                                    
param1_2 = 1E8; %Peak Frequency Deviaton in Hz
param2_2 = 1E5; %Modulation frequency in HZ

% Parameters for Digital Modulation
numOfSymbols2 = 2048;
symbolRate2 = 30E6;
rollOff2 = 0.35;

% IQ Complex Waveform Calculation
if modType1 <= 4
    wfmIq = CalculateAnalogModWfm(  modType1,...
                                    minCycles,...
                                    samplingRate,...
                                    interpol / intCorrection,...
                                    granul / gCorrection,...
                                    param1_1,...
                                    param2_1);
else
    wfmIq = CalculateDigitalModWfm( modType1,...
                                    numOfSymbols1,...                        
                                    symbolRate1,...
                                    rollOff1,... 
                                    samplingRate,... 
                                    interpol / intCorrection);   
    % wfmIq length is not adjusted for granularity to optimize accuracy for
    % symbol rate. Howeverm symbol rate must be adjusted for signal loop
    % consistency. Actual Symbol Rate must be calculated and used in 
    % analysis of the signal.                          
    actualSymbR = samplingRate / interpol * numOfSymbols1 / length(wfmIq);
    fprintf('\nActual Symbol Rate for DUC1 = %d\n', actualSymbR);
    
end

% I and Q waveforms
myWfmI1 = real(wfmIq);
myWfmQ1 = imag(wfmIq);
clear wfmIq;
% Negative Q wavefomr for inverse spectrum
if reverse
    myWfmQ1 = -myWfmQ1;
end

% Frequency Offset applied
[myWfmI1,  myWfmQ1] = ApplyFreqOffset(  fOffset1,...
                                        samplingRate / interpol,...
                                        myWfmI1,...
                                        myWfmQ1);
                                    


% IQ Complex Waveform Calculation for second DUC
if modType2 <= 4
    wfmIq = CalculateAnalogModWfm(  modType2,...
                                    minCycles,...
                                    samplingRate,...
                                    interpol / intCorrection,...
                                    granul / gCorrection,...
                                    param1_2,...
                                    param2_2);
else
    wfmIq = CalculateDigitalModWfm( modType2,...
                                    numOfSymbols2,...                        
                                    symbolRate2,...
                                    rollOff2,... 
                                    samplingRate,... 
                                    interpol / intCorrection);
                                
    actualSymbR = samplingRate / interpol * numOfSymbols2 / length(wfmIq);
    fprintf('\nActual Symbol Rate for DUC2 = %d\n', actualSymbR);
    
end

myWfmI2 = real(wfmIq);
myWfmQ2 = imag(wfmIq);
clear wfmIq;
% Negative Q wavefomr for inverse spectrum
if reverse
    myWfmQ2 = -myWfmQ2;
end

% Frequency Offset applied
[myWfmI2,  myWfmQ2] = ApplyFreqOffset(  fOffset1,...
                                        samplingRate / interpol,...
                                        myWfmI2,...
                                        myWfmQ2);
% Joint normalization for each IQ pair so clipping is avoided                                    
[myWfmI1,  myWfmQ1] = NormalIq(myWfmI1, myWfmQ1);
[myWfmI2,  myWfmQ2] = NormalIq(myWfmI2, myWfmQ2);
% Waveform data formated in a single uint8 array for downooad as a single
% segment.
myWfm = formatWfm2(myWfmI1, myWfmQ1, myWfmI2, myWfmQ2);
clear myWfmI1 myWfmI2;
clear myWfmQ1 myWfmQ2;

% If necessary, wfm repetated for waveform granularity
myWfm = trimGranTwo(myWfm, granul / 2);

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Reset AWG
inst.SendCmd('*CLS');
inst.SendCmd('*RST');

% Set sampling rate for AWG to maximum.
inst.SendCmd([':FREQ:RAST ' num2str(2.5E9)]);
inst.SendCmd(sprintf(':INST:CHAN %d', channel));
% Interpolation factor for I/Q waveforms set to X8 but it is x16
inst.SendCmd(':SOUR:INT X8');
% DAC Mode set to 'DUC' and IQ Modulation mode set to 'TWO'
inst.SendCmd(':MODE DUC');
inst.SendCmd(':IQM TWO'); 

% Waveform Downloading
% *******************
inst.SendCmd(':TRAC:DEL:ALL');
fprintf(1, 'DOWNLOADING WAVEFORM\n');
res = SendWfmToProteusTwo(inst, channel, segment, myWfm);
fprintf(1, 'WAVEFORM DOWNLOADED!\n');
clear myWfm;

% Select segment for generation
fprintf(1, 'SETTING AWG OUTPUT\n');
inst.SendCmd(sprintf(':FUNC:MODE:SEGM %d', segment))
% Output volatge set to MAX
inst.SendCmd(':SOUR:VOLT MAX');
% 6dB IQ Modulation gain applied
if apply6db
    inst.SendCmd(':NCO:SIXD2 ON');   
else
    inst.SendCmd(':NCO:SIXD2 OFF');    
end
% NCO frequency and phase setting
inst.SendCmd(sprintf(':NCO:CFR1 %d', ncoFreq2));
inst.SendCmd(sprintf(':NCO:PHAS1 %d', phase2));
% 6dB IQ Modulation gain applied
if apply6db
    inst.SendCmd(':NCO:SIXD1 ON');   
else
    inst.SendCmd(':NCO:SIXD1 OFF');    
end
% NCO frequency and phase setting
inst.SendCmd(sprintf(':NCO:CFR2 %d', ncoFreq1));
inst.SendCmd(sprintf(':NCO:PHAS2 %d', phase1));
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

function finalWfm = trimGranTwo(inWfm, granularity)
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
    
    finalWfm = uint8(zeros(1, finaL));
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

function [normI,  normQ] = NormalIq(wfmI, wfmQ)    
    
    maxPwr = max(wfmI.*wfmI + wfmQ .* wfmQ);
    maxPwr = maxPwr ^ 0.5;
    
    normI = wfmI / maxPwr;
    normQ = wfmQ / maxPwr;
end

function outWfm = Interleave2(wfmI, wfmQ)   

    wfmLength = length(wfmI);
    if length(wfmQ) < wfmLength
        wfmLength =  length(wfmQ);
    end
    
    %wfmLength = 2 * wfmLength;
    outWfm = uint8(zeros(1, 2 * wfmLength));
    
    outWfm(1:2:(2 * wfmLength - 1)) = wfmI;
    outWfm(2:2:(2 * wfmLength)) = wfmQ;
end

function retval = myQuantization (myArray, dacRes)
  
  minLevel = 1; % To avoid DC component
  maxLevel = 2 ^ dacRes - 1;  
  numOfLevels = maxLevel - minLevel + 1;
  
  retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
  retval = retval + minLevel;
  
  retval(retval > maxLevel) = maxLevel;
  retval(retval < minLevel) = minLevel;

end

function outWfm = formatWfm2(inWfmI1, inWfmQ1, inWfmI2, inWfmQ2)
%formatWfm2 This function formats data for two I/Q streams to be dwnloaded
%to a single segment in Proteus to be generated in the IQM Mode 'TWO'
%   All waveforms must be properly normalized to the -1.0/+1.0 range.
%   All waveforms must have the same length

    % Formatting requires to go through the following steps:
    %   1) quantize samples to 16-bit unsigned integers
    %   2) swap the LSB and MSB as MSBs will be sent first for this mode
    %   3) convert the uint16 array to an uint8 array of twice the size
    % Final wfm is MSB, LSB, MSB, LSB,...
    inWfmI1 = typecast(swapbytes(uint16(myQuantization(inWfmI1, 16))),'uint8'); 
    inWfmQ1 = typecast(swapbytes(uint16(myQuantization(inWfmQ1, 16))),'uint8');
    inWfmI2 = typecast(swapbytes(uint16(myQuantization(inWfmI2, 16))),'uint8');
    inWfmQ2 = typecast(swapbytes(uint16(myQuantization(inWfmQ2, 16))),'uint8');
    % Sequence MSBI1, MSBQ1, MSBQ2, MSBI2, LSBI1, LSBQ1, LSBQ2, LSBI2
    % This is done in three interleaving steps
    outWfmI = Interleave2(inWfmI1, inWfmQ2);
    outWfm = Interleave2(inWfmQ1, inWfmI2);
    outWfm = Interleave2(outWfmI, outWfm);    
end

function result = SendWfmToProteusTwo(  instHandle,...
                                        channel,...
                                        segment,...
                                        myWfm)

    %Select Channel
    instHandle.SendCmd(sprintf(':INST:CHAN %d', channel));  
    % Space rquired by the segment is half the length of the uinit8 input
    % array as segments are defined as 16-bit samples
    instHandle.SendCmd(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)/2));        
    % select segmen as the the programmable segment
    instHandle.SendCmd(sprintf(':TRAC:SEL %d', segment));
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';    
    instHandle.SendBinaryData(prefix, myWfm, 'uint8');       
    
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

    