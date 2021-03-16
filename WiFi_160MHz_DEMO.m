% 802.11ax DEMO including Envelope Tracking
%% Generate 802.11ax Waveform

%802.11ax Definition
chBw = 'CBW160';
apepLength = 10000;
numPackets = 1;
idleTime = 100E-6;

%AWG Generation Settings
samplingRate = 8.96E9;
ampl = 0.8;
interpol = 8;
dataSize = 16;
granul = 64;
reverse = false;

cfr = (2.412) * 1E+09;%2.412e9;
fOffset = 0.0E6;
phase = 0.0 ;
channel = 1;
segment = 1;
should_reset = true;
apply6db = true;

% Envelope Tracking
genEnvlp = true;
minLvl = 0.1;


% 802.11ax configuration:
heSUCfg = wlanHESUConfig('ChannelBandwidth', chBw, ...
    'NumTransmitAntennas', 1, ...
    'NumSpaceTimeStreams', 1, ...
    'SpatialMapping', 'Direct', ...
    'PreHESpatialMapping', false, ...
    'MCS', 5, ...
    'ChannelCoding', 'LDPC', ...
    'APEPLength', apepLength, ...
    'GuardInterval', 3.2, ...
    'HELTFType', 4, ...
    'UplinkIndication', false, ...
    'BSSColor', 0, ...
    'SpatialReuse', 0, ...
    'TXOPDuration', 127, ...
    'HighDoppler', false);


% input bit source:
in = randi([0, 1], 1000, 1);


% waveform generation:
waveform = wlanWaveformGenerator(in, heSUCfg, ...
    'NumPackets', numPackets, ...
    'IdleTime', idleTime, ...
    'ScramblerInitialization', 93, ...
    'WindowTransitionTime', 1e-07);

Fs = wlanSampleRate(heSUCfg); 								 % sample rate of waveform

%% Visualize 802.11ax Waveform
% Time Scope
%timeScope = dsp.TimeScope('SampleRate', Fs, ...
%    'TimeSpanOverrunAction', 'Scroll', ...
%    'TimeSpan', 1.875e-07);
%timeScope(waveform);
%release(timeScope);

% Spectrum Analyzer
%spectrum = dsp.SpectrumAnalyzer('SampleRate', Fs);
%spectrum(waveform);
%release(spectrum);

% RU Assignment and Allocated Subcarriers:
%showAllocation(heSUCfg);

% Rate Conversion and formatting

dataI = real(waveform);
dataQ = imag(waveform);
clear wfmIq;

dataI = dataI';
dataQ = dataQ';

% waveform generation:
% Actual waveform sample rate calculation
samplingRateBb = samplingRate / interpol;
% Fractional reduction of sample rate versus symbol rate
% for fractional resampling
[num, den] = reduceFraction(samplingRateBb, Fs);
fprintf(1, 'RESAMPLING BASEBAND WAVEFORMS\n');
finalI = MyIdealInterpolation2(dataI, num, 100000, 1);
finalQ = MyIdealInterpolation2(dataQ, num, 100000, 1);

clear dataI dataQ;

% The basic waveform must be repeated until the total number of samples
% is a multiple of the decimation factor (den)to kepp consistency
if den > 1
    fprintf(1, 'DECIMATING BASEBAND\n');
    finalI = trimGran(finalI, den);
    finalQ = trimGran(finalQ, den);
    % Actual decimation: one every den samples is actually preserved
    finalI = finalI(1:den:length(finalI));
    finalQ = finalQ(1:den:length(finalQ));
end

fprintf(1, 'NORMALIZING BASEBAND WAVEFORMS\n');
[finalI, finalQ, divFactor] = myNormalization3 (finalI, finalQ, true);

fprintf(1, 'ADJUSTING GRANULARITY\n');
finalI = trimGran(finalI, granul);
finalQ = trimGran(finalQ, granul);

if genEnvlp
    %waveform = waveform';
    envlpWfm = GetEnvlpWfm(finalI + 1i * finalQ, minLvl);
end

% Final signal conditioning

if reverse
    finalQ = -finalQ;
end

[finalI,  finalQ] = ApplyFreqOffset(fOffset, samplingRate / 8, finalI, finalQ);

fprintf(1, 'I/Q INTERLEAVING\n');
%[finalI,  finalQ] = NormalIq(finalI, finalQ); 
myWfm = Interleave(finalI, finalQ);
clear finalI finalQ;

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 
fprintf('Connecting to Proteus\n');
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

%% Create Administrator
% Create instance of the CProteusAdmin class, and open it.
% Note that only a single CProteusAdmin instance can be open, 
% and it must be kept open till the end of the session.

admin = CProteusAdmin(@OnLoggerEvent);

rc = admin.Open();
assert(rc == 0);


slotIds = admin.GetSlotIds();
numSlots = length(slotIds);
assert(numSlots > 0);
% If there are multiple slots, let the user select one ..
sId = slotIds(1);
if numSlots > 1
    fprintf(1, '\n%d slots were found\n', numSlots);
    for n = 1:numSlots
        sId = slotIds(n);
        slotInfo = admin.GetSlotInfo(sId);
        if ~slotInfo.IsSlotInUse
            modelName = slotInfo.ModelName;
            if slotInfo.IsDummySlot
                fprintf(1, ' * Slot Number: Model %s [Dummy Slot].\n', sId, ModelName);
            else
                fprintf(1, ' * Slot Number: Model %s.\n', sId, ModelName);             
            end
        end
    end
    choice = input('Enter SlotId ');
    fprintf(1, '\n');
    sId = uint32(choice);
end

% Connect to the selected instrument ..
inst = admin.OpenInstrument(sId, should_reset);
instId = inst.InstrId;

res = inst.SendScpi('*IDN?');
assert(res.ErrCode == 0);
fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));

awgId = netStrToStr(res.RespStr);

fprintf(1, 'SETTING AWG\n');

res = inst.SendScpi('*CLS');
assert(res.ErrCode == 0);

res = inst.SendScpi('*RST');
assert(res.ErrCode == 0);

% Set sampling rate for AWG as defined in the preamble.
res = inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));


% ---------------------------------------------------------------------
% Play selected segment in selected channels and set amplitude
% ---------------------------------------------------------------------

res = inst.SendScpi(sprintf(':INST:CHAN %d', channel));
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

res = inst.SendScpi(':TRAC:DEL:ALL');
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

res = inst.SendScpi(':SOUR:MODE DUC');
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

res = inst.SendScpi(':SOUR:INT X8');
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

res = inst.SendScpi(':SOUR:IQM ONE'); 
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

fprintf(1, 'DOWNLOADING WAVEFORM\n');
%*******************************************************
res = SendWfmToProteus(inst, channel, segment, myWfm, 16);
assert(res.ErrCode == 0);
fprintf(1, 'WAVEFORM DOWNLOADED!\n');
clear myWfm;

res = inst.SendScpi(':TRAC:FORM?');
assert(res.ErrCode == 0);
fprintf(1, '\nTRACE FORMAT: ''%s''\n', netStrToStr(res.RespStr));

fprintf(1, 'SETTING AWG OUTPUT\n');
res = inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segment));
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

res = inst.SendScpi(sprintf(':SOUR:CFR %d', cfr));
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

res = inst.SendScpi(sprintf(':SOUR:PHAS %d', phase));
assert(res.ErrCode == 0);

res = inst.SendScpi(sprintf(':SOUR:VOLT %d', ampl));
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

if apply6db
    res = inst.SendScpi(':SOUR:SIXD ON');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
else
    res = inst.SendScpi(':SOUR:SIXD OFF');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
end

res = inst.SendScpi(':OUTP ON');
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
fprintf(1, 'SETTING SAMPLING CLOCK\n');
% Set sampling rate for AWG as defined in the preamble.
res = inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);
assert(res.ErrCode == 0);

res = inst.SendScpi(':SYST:ERR?');
assert(res.ErrCode == 0);
fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

%Envelope Tracking
%*****************************************

if genEnvlp

    envlpWfm = Interleave(envlpWfm, zeros(1, length(envlpWfm)));

    % ---------------------------------------------------------------------
    % Play selected segment in selected channels and set amplitude
    % ---------------------------------------------------------------------

    res = inst.SendScpi(sprintf(':INST:CHAN %d', channel + 1));
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    res = inst.SendScpi(':SOUR:MODE DUC');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    res = inst.SendScpi(':SOUR:INT X8');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    res = inst.SendScpi(':SOUR:IQM ONE'); 
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    fprintf(1, 'DOWNLOADING ENVELOPE WAVEFORM\n');
    %*******************************************************
    res = SendWfmToProteus(inst, channel + 1, segment + 1, envlpWfm, 16);
    assert(res.ErrCode == 0);
    fprintf(1, 'WAVEFORM DOWNLOADED!\n');
    clear myWfm;

    res = inst.SendScpi(':TRAC:FORM?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nTRACE FORMAT: ''%s''\n', netStrToStr(res.RespStr));

    fprintf(1, 'SETTING AWG OUTPUT\n');
    res = inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segment + 1));
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    res = inst.SendScpi(sprintf(':SOUR:CFR %d', 0.0));
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    res = inst.SendScpi(sprintf(':SOUR:PHAS %d', 0.0));
    assert(res.ErrCode == 0);

    res = inst.SendScpi(sprintf(':SOUR:VOLT %d', ampl));
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));

    if apply6db
        res = inst.SendScpi(':SOUR:SIXD ON');
        assert(res.ErrCode == 0);

        res = inst.SendScpi(':SYST:ERR?');
        assert(res.ErrCode == 0);
        fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    else
        res = inst.SendScpi(':SOUR:SIXD OFF');
        assert(res.ErrCode == 0);

        res = inst.SendScpi(':SYST:ERR?');
        assert(res.ErrCode == 0);
        fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    end

    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    fprintf(1, 'SETTING SAMPLING CLOCK\n');
    % Set sampling rate for AWG as defined in the preamble.
    res = inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
end




%**************************************

% It is recommended to disconnect from instrumet at the end
rc = admin.CloseInstrument(instId);   
% Close the administrator at the end
admin.Close();

clear;
fprintf(1, 'END\n');

function envlpWfm = GetEnvlpWfm(inWfm, minLvl)
    envlpWfm = abs(inWfm);
    envlpWfm(envlpWfm < minLvl) = minLvl;
end

function finalWfm = trimGran(inWfm, granularity)

    baseL = length(inWfm);
    finaL = lcm(baseL, granularity);
    
    finalWfm = zeros(1, finaL);
    pointer = 1;
    
    while pointer < finaL
        finalWfm(pointer : (pointer + baseL -1)) = inWfm;
        pointer = pointer + baseL;        
    end

end

function retval = MyIdealInterpolation2 (myArray, xFactor, quality, bwFraction)
  
  %expansion by zero-padding
  retval = zeros(1, xFactor * length(myArray));
  retval([1:xFactor:end]) = myArray;
  % "Ideal" Interpolation filter
  lenSinc = quality; 
  
  mySinc = -lenSinc * xFactor : 1 : lenSinc * xFactor ;  
  mySinc = sinc(bwFraction .* mySinc / xFactor);
  %The FIR taps must be normalized so the DC response is 0dB
  normFactor = sum(mySinc(1:xFactor:length(mySinc)));
  normFactor = 1/ normFactor;
  mySinc = normFactor .* mySinc;
  myWindow = blackman(length(mySinc));
  myWindow = myWindow.'; 
  mySinc = mySinc .* myWindow;
  %convolution
  retval = cconv(retval, mySinc, length(retval));
  %retval = real(retval);
  
end

function [normI, normQ, divFactor] = myNormalization3 (iArray, qArray, keepZero)
  
  maxLevel = max(iArray);
  if max(qArray) > maxLevel
      maxLevel = max(qArray);
  end
  
  minLevel = min(iArray);
  if min(qArray) < minLevel
      minLevel = min(qArray);
  end
  
  a = 1; 
  b = 0;
  
  if keepZero      
      if abs(maxLevel) > 0 || abs(minLevel) > 0
          a = abs(maxLevel);         
          if (abs(minLevel)> a)
              a = abs(minLevel);  
          end      
          a = 1 / a;         
      end
  else
      if maxLevel > minLevel
          a = 0.5 * (maxLevel - minLevel);
          a = 1 / a;
          b = maxLevel + minLevel;
          b = b / 2;    
      end
  end
  
  a = 1.0 * a; %Reduce the 1.0 if intrpolation filter causes clipping
  normI = a * (iArray - b);
  normQ = a * (qArray - b);
  divFactor = a;

end

function [outNum, outDen] = reduceFraction(num, den)
%reduceFraction Reduce num/den fraction
%   Use integers although not mandatory
    num = round(num);
    den = round(den);
    % Reduction is obtained by calcultaing the greater common divider...
    G = gcd(num, den);
    % ... and then dividing num and den by it.
    outNum = num / G;
    outDen = den / G;
end

% Function netStrToStr
function str = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch        
        str = '';
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

function outWfm = Interleave(wfmI, wfmQ)    
    wfmLength = length(wfmI);
    if length(wfmQ) < wfmLength
        wfmLength =  length(wfmQ);
    end
    
    %wfmLength = 2 * wfmLength;
    outWfm = zeros(1, 2 * wfmLength);
    
    outWfm(1:2:(2 * wfmLength - 1)) = wfmI;
    outWfm(2:2:(2 * wfmLength)) = wfmQ;
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

function result = SendWfmToProteus(instHandle, channel, segment, myWfm, dacRes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
    %Select Channel
    res = instHandle.SendScpi(sprintf(':INST:CHAN %d', channel));
    assert(res.ErrCode == 0);   
    
    res = instHandle.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    
    % Set sampling rate for AWG as defined in the preamble.
    if dacRes == 8
        res = instHandle.SendScpi(':TRAC:FORM U8');
        assert(res.ErrCode == 0);
    else
        res = instHandle.SendScpi(':TRAC:FORM U16');
        assert(res.ErrCode == 0);
    end

    res = instHandle.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    
    res = instHandle.SendScpi(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)));
    assert(res.ErrCode == 0);
    
    res = instHandle.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    
    % select segmen as the the programmable segment
    res = instHandle.SendScpi(sprintf(':TRAC:SEL %d', segment));
    assert(res.ErrCode == 0); 
    
    res = instHandle.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    
    % format Wfm
    myWfm = formatWfm(myWfm,dacRes);
    
    % Download the binary data to segment    
    res = instHandle.WriteBinaryData([':TRAC:DATA 0,#' num2str(length(myWfm))], myWfm);
    assert(res.ErrCode == 0);
    
    res = instHandle.SendScpi('*OPC?');
    assert(res.ErrCode == 0);
    
    res = instHandle.SendScpi(':SYST:ERR?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nCommand Error: ''%s''\n', netStrToStr(res.RespStr));
    
    result = res;
end

