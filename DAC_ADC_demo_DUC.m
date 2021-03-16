% DAC/ADC Proteus Demo with DUC
%clear all;
close all;
clear variables;
clear global;
clc;
paranoia_level = 0; % 0, 1 or 2
fprintf('ADC/DAC DEMO STARTS\n');
% Signal Parameters
c   = 3e8;
fprintf('\nParameter Definition\n');
PW  = [10.0E-6 10.0E-6 10.0E-6 10.0E-6 10.0E-6 10.0E-6 10.0E-6];
BW  = [5.0E+06 5.0E+06 5.0E+06 5.0E+06 5.0E+06 5.0E+06 5.0E+06];
PRF = [4.0E+03 8.0E+03 16.0E+03 4.0E+03 4.0E+03 4.0E+03 4.0E+03];
FC  = [3.2E+09 3.2E+09 3.2E+09 3.2E+09 1.2E+09 0.5E+09 3.2E+09];
NP  = [16 32 64 16 16 16 16];
DL  = [0.0 0.0 0.0 0.0 0.0 0.0 100.0E-06];
DS  = [0.0 0.0 0.0 0.0 0.0 0.0 100];

% Signal to be generated
NS  = 1;
numOfAcq = 11;

numberOfFrames = 1;

cFreq = FC(NS);
span = BW(NS);

% AWG Parameters

connStr = '192.168.1.48'; % your IP here

samplingRate    = 9E+09; %Sampling rate for the target AWG
dacRes          = 16; %Default DAC resolution
granularity     = 64; %Waveform length granularity
channelNumber   = 1;  % Channel number to generate mmlti-tone
segmentNumber   = 1;  %Segment # to store the multi-tone wfm
voltLvl         = 0.5; % Ampltiude level for the AWG
interpol        = 8;
% Digitizer Parameters

samplingRateDig = 2.4E+09;
adcGranul       = 48;
dacResDig       = 12;
preTrig = 4800;
frameLen = samplingRateDig * PW(NS);
frameLen = frameLen + 2 * preTrig;
frameLen = ceil(frameLen / adcGranul) * adcGranul;

acqWfm = zeros(numOfAcq, frameLen); 


cFreq = FC(NS);
while cFreq > samplingRateDig / 2
    cFreq = cFreq - samplingRateDig;
    if cFreq < 0
        cFreq = -cFreq;
    end
end
span = 4 * BW(NS);

% Connecting to Proteus
%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 
fprintf('Connecting to Proteus\n');
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

%% Create Administrator and Select Proteus Unit
% Create instance of the CProteusAdmin class, and open it.
% Note that only a single CProteusAdmin instance can be open, 
% and it must be kept open till the end of the session.

try
    %connStr = input('Enter the instrument''s IP address  ', 's');
    
    inst = TEProteusInst(connStr, paranoia_level);
    fprintf('\n');
    
    res = inst.Connect();
    assert (res == true);
    
    idnstr = inst.SendQuery('*IDN?');
    fprintf('Connected to: %s\n', idnstr);
    
    fprintf('Reset instrument ..\n');
    inst.SendCmd('*CLS; *RST');
    
    % Set sampling rate for AWG as defined in the preamble.   
    inst.SendCmd([':FREQ:RAST ' num2str(samplingRate)]);    
    
    inst.SendCmd(['INST:CHAN ' num2str(channelNumber)]);
    
    inst.SendCmd(':TRAC:DEL:ALL');
   
    inst.SendCmd(':SOUR:MODE DUC');    

    inst.SendCmd(':SOUR:INT X8');

    inst.SendCmd(':SOUR:IQM ONE'); 
    
    % Signal Creation Starts Here
    % ***************************
    actualSR = samplingRate / interpol;
    fprintf('\nChirp Waveform Calculation\n');
    totDuration = NP(NS) / PRF(NS) + DL(NS);
    totDuration = totDuration * c / (c + DS(NS)); %doppler
    numOfSamples = totDuration * actualSR;
    numOfSamples = granularity * ceil(numOfSamples / granularity);

    myWfm = zeros(1, numOfSamples);
    mywfm = complex(myWfm);

    %isolated I/Q complex pulse
    isoPulse = r_wf(PW(NS), BW(NS), actualSR, DS(NS), 1);
    lPulse = length(isoPulse);

    % baseband isolated pulses are added to the wfm
    for i = 0:(NP(NS) - 1)
        pulseLocation = DL(NS) + i / PRF(NS);
        pulseLocation = pulseLocation * c / (c + DS(NS)); % Doppler
        pulseLocation = round(pulseLocation * actualSR);
        myWfm((pulseLocation + 1):(pulseLocation + lPulse)) = isoPulse;    
    end
    %plot(imag(myWfm));
    % carrier is modulated
    dopplerF = FC(NS) * (c + DS(NS))/c - FC(NS);
    
    
    myWfm = ApplyFreqOffset2(dopplerF, samplingRate / 8, myWfm);
    %myWfm = IqModulator2 (myWfm, dopplerF, 0, actualSR); 
    
    myWfm = Interleave(real(myWfm),imag(myWfm));
    
    %hold;  

    % Waveform will be downoaded later.
    fprintf('\nChirp Waveform Download to AWG\n');   

    myWfmLength = length(myWfm);
    res = SendWfmToProteus(inst, channelNumber, segmentNumber, myWfm, 16);     

    fprintf('Chirp Waveform Downloaded\n');
    
    % ---------------------------------------------------------------------
    % Play selected segment in selected channel and set amplitude
    % ---------------------------------------------------------------------

    inst.SendCmd(sprintf(':SOUR:CFR %d', FC(NS)));
    
    inst.SendCmd(':SOUR:SIXD ON');
    
    inst.SendCmd([':SOUR:FUNC:MODE:SEGM ' num2str(segmentNumber)]);
   

    inst.SendCmd([':SOUR:VOLT ' num2str(voltLvl)]);
 

    inst.SendCmd(':OUTP ON');
    
    % Set sampling rate for AWG as defined in the preamble.   
    inst.SendCmd([':FREQ:RAST ' num2str(samplingRate)]);

    % -------------------------------------
    % Operate the ADC with direct functions
    % -------------------------------------    
       
    ddcFactor = 32;
    
    matchedFilter = GetMatchedFilter(FC(NS), PW(NS), BW(NS), samplingRateDig);
    matchedFilter = matchedFilter(1:ddcFactor:length(matchedFilter));
    
    % ----------
    % ADC Config
    % ----------
    
    % Turn on ADC dual-channels mode
    inst.SendCmd(':DIG:MODE DUAL');

    % Set sampling rate
    
    inst.SendCmd(':DIG:FREQ %g', samplingRateDig);
    
    for chan = 1:1
        % Select digitizer channel:
        inst.SendCmd(':DIG:CHAN %d', chan);
        % Set the voltage-range of the selected channel to maximum (1V)
        inst.SendCmd(':DIG:CHAN:RANG MED');
    end
    
    % Enable acquisition in the digitizer's channels
    %inst.SendCmd(':DIG:CHAN 2; :DIG:CHAN:STATE ENAB');
    inst.SendCmd(':DIG:CHAN 1; :DIG:CHAN:STATE ENAB');
    
    % Setup frames layout
    
    inst.SendCmd(':DIG:ACQ:DEF %d, %d',numberOfFrames, frameLen);
    
    % Set channel 1 of the digitizer as it trigger source
    inst.SendCmd(':DIG:TRIG:SOUR CH1');
    inst.SendCmd(':DIG:PRET %d', preTrig);

    % Set the threshold level of the self-trigger (of channel 1):

    triglevel = 0.025;
   
    % Select channel 1
    inst.SendCmd(':DIG:CHAN 1');

    % Set the threshold level
    inst.SendCmd(':DIG:TRIG:SELF %f', triglevel);
    
    
    % Select which frames are filled with captured data 
    %(all frames in this example)
    inst.SendCmd(':DIG:ACQ:CAPT:ALL');
    

    for i=1:numOfAcq
        if i== 2
            tic
        end
        
        inst.SendCmd(':DIG:INIT ON');    
        
        for n = 1:250
            resp = inst.SendQuery(':DIG:ACQ:STATus?');
            resp = strtrim(resp);
            items = split(resp, ',');
            items = str2double(items);
            if length(items) >= 3 && items(2) == 1
                break
            end
            if mod(n, 10) == 0
                fprintf('%d. %s\n', fix(n / 10), resp);
            end
            pause(0.1);
        end
       
        %resp = inst.SendQuery(':DIG:ACQ:STATus?');
        %resp = strtrim(resp);

        inst.SendCmd(':DIG:INIT OFF');

        % Define what we want to read 
        % (frames data, frame-header, or both).
        % In this example we read the frames-data
        inst.SendCmd(':DIG:DATA:TYPE FRAM');
        inst.SendCmd(':DIG:DATA:SEL ALL');

        %fprintf('Total size of %d frames in bytes %d\n', numberOfFrames, totSize);

        % Read binary block
        samples = inst.ReadBinaryData(':DIG:DATA:READ?', 'uint16');
        samples = int16(samples) - 2048;
        acqWfm(i,:) = samples;
        samples2 = samples(1:ddcFactor:length(samples));
        
        % Pulse Compression
  
        myWfm = r_conv2(samples2, matchedFilter);
        %fprintf('Pulse Compression Completed\n');
        
    end
    
    toc
    
    itemToShow = numOfAcq;
    ShowResults(acqWfm, matchedFilter, itemToShow, preTrig, cFreq, span, ddcFactor, samplingRateDig);
    % ---------------------------------------------------------------------
    % End of the example
    % ---------------------------------------------------------------------
       
    % Close the session
    inst.Disconnect();
    
    clear inst
  
    fprintf('\nELBIT DEMO COMPLETED\n');
    
    fig = uifigure('Position',[100 100 350 100]);    

    sld = uislider(fig,...
        'Position',[10 75 300 3],...
        'Limits', [1 numOfAcq],...
        'ValueChangedFcn',@(sld,event) ShowResults( acqWfm,...
                                                    matchedFilter,...
                                                    sld.Value,...
                                                    preTrig,...
                                                    cFreq,...
                                                    span,...
                                                    ddcFactor,...
                                                    samplingRateDig));

catch ME    
    close all;
    rethrow(ME)
end

function ShowResults(acqWfm, matchedFilter, itemToShow, preTrig, cFreq, span, ddcFactor, samplingRateDig)

    itemToShow = round(itemToShow);
    
    if itemToShow < 1
        itemToShow = 1;
    end
    
    if itemToShow > length(acqWfm(:,1))
        itemToShow = length(acqWfm(:,1));
    end
    
    samples = acqWfm(itemToShow,:);
    samples2 = samples(1:ddcFactor:length(samples));
    myWfm = r_conv2(samples2, matchedFilter);
    
    tiledlayout(2,2);
    
    % Top plot
    ax1 = nexttile;
    xData = 0:(length(samples) - 1);
    xData = xData - preTrig;
    xData = xData / samplingRateDig;
    
    xData = xData * 1e+06;
    
    plot(ax1, xData, samples)
    title(ax1,sprintf('Acquired Waveform # %d', itemToShow));
    xlabel('microseconds') 
    ylabel('ADC Levels') 
    grid(ax1,'on')
    
    % Mid plot
    ax3 = nexttile;
    pSpec = abs(fft(samples));
    
    
    startF = cFreq - span/2;
    stopF = startF + span;
    
    xData = 0:(length(samples) - 1);
    xData = xData * samplingRateDig / length(samples);
    
    xData1 = xData(find(xData >= startF & xData <= stopF));
    pSpec = pSpec(find(xData >= startF & xData <= stopF));
    pSpec = 20 * log10(pSpec / max(pSpec));
    
    xData1 = xData1 / 1e+06;
    
    plot(ax3, xData1, pSpec);
    title(ax3,'Acquired Waveform Spectrum')
    xlabel('MHz') 
    ylabel('dB') 
    grid(ax3,'on')
    %Plot Pulse Compression signal 
    % Bottom plot
    nexttile([1 2]);
    %ax2 = nexttile;
    xData = 0:(length(myWfm) - 1);
    xData = ddcFactor * xData / samplingRateDig;
    
    xData = xData * 1e+06;
    
    plot(xData, myWfm, 'LineWidth',1)
    title('Pulse Compression')
    xlabel('microseconds') 
    ylabel('Level') 
    grid('on')


end

% Function netStrToStr
function str = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch        
        str = '';
    end
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

function outWfm = formatWfm(inWfm,dacRes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    outWfm = myQuantization(inWfm, dacRes); 
    if dacRes == 8
        outWfm = uint8(outWfm);
    elseif dacRes == 16
        outWfm = uint16(outWfm);
        outWfm = typecast(outWfm, 'uint8');
    end
end

function retval = myNormalization (myArray, keepZero)
  
  maxLevel = max(myArray);
  minLevel = min(myArray);
  
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
  
  retval = a * (myArray - b);

end

function result = SendWfmToProteus(instHandle, channel, segment, myWfm, dacRes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
    %Select Channel
    instHandle.SendCmd(sprintf(':INST:CHAN %d', channel));
   
    
    instHandle.SendCmd(':SYST:ERR?');
   
    instHandle.SendCmd(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)));
        
    % select segmen as the the programmable segment
    instHandle.SendCmd(sprintf(':TRAC:SEL %d', segment));

    % format Wfm
    %myWfm = formatWfm(myWfm,dacRes);
    myWfm = myQuantization(myWfm, 16);
    
    % Download the binary data to segment  
    
    prefix = ':TRAC:DATA 0';
    
    instHandle.SendBinaryData(prefix, myWfm, 'uint16');
    
    result = length(myWfm);
end



function samples = GetAdcWfm(inst, ch, sampleRate, frameLen)
%UNTITLED3 Summary of this function goes here
    
    numFrames = 1;
    offLen = 0;
    
    rc = inst.SetAdcAcquisitionEn(1,1);
    assert(rc == 0);
    
    rc = inst.SetAdcFramesLayout(numFrames, frameLen);
    assert(rc == 0);  

    % ---------------------------------------------------------------------
    % Capture Parmeters
    % ---------------------------------------------------------------------
       
    rc = inst.SetAdcMultiFrameModeEn(1);
    assert(rc == 0);
    % ---------------------------------------------------------------------
    % ADC Capture
    % ---------------------------------------------------------------------
              
    rc = inst.SetAdcCaptureEnable(1);
    assert(rc == 0);
    %choice = input('Wait for shuttle to complete then press [Enter] to continue ');
    %fprintf(1, '\n');
    % Wait till the capture completes
    
    status = inst.ReadAdcCaptureStatus();
    for i = 1 : 250
        if status ~= 0
            break;
        end
        pause(0.01);
        status = inst.ReadAdcCaptureStatus();
    end
    
    rc = inst.SetAdcCaptureEnable(0);
    assert(rc == 0);
    
    %chanIndex = 0;
	frame1st = 1;
    numFrames = 1;
    adcChanInd = ch - 1;
    
	% allocate NET array
    netArray = NET.createArray('System.UInt16', frameLen * numFrames);
    
    % read the captured frame
    rc = inst.ReadMultipleAdcFrames(adcChanInd, frame1st, numFrames, netArray);

	% cast to matlab vector
	samples = double(netArray);
    
	% deallocate the NET array
    delete(netArray);
    
    figure(1); clf;
    samples = samples - 2048;
    %plot(samples);
    samples = myNormalization(samples, 1);
    %plot(samples);      
    % Free the memory space that was allocated for ADC capture
    rc = inst.FreeAdcReservedSpace();
    assert(rc == 0);
end

function x = r_wf(tau, Beta, Fs, dopplerShift, wf_type)
    % This function :Create waveform to Signal Generator
    % 1. Creates a LFM (IQ) or  chirp (IF)  
    % 2. Plots signal and records to file
    % 3. Specs: tau=pulse width, Beta=frequency sweep range, % Fs = sampling rate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Ts      = 1/Fs;
    c       = 3e8;
    tau     = tau * c / (c + dopplerShift);
    N       = ceil(tau*Fs);  %Number of sample points of input chirp signal
    t       = 0:Ts:Ts*(N-1);
    t0      = abs(tau/2);    
    fc      = 0;%100.0e6; %(base band freq)

    %%sig_pulse=A*exp(1j*2*pi*(160e6-fc)*T).*sig;

    if wf_type==1
    % LFM complex IQ signal    
        mu  = Beta/tau;
        corrTime = (t-t0) * (c + dopplerShift) / c;
        x   = exp(1j*pi*mu*(corrTime.^2));    
    else
    % Chirp IF         
        x  = chirp(t,fc,tau,Beta);    
    end     
end

function [filter_out]=r_conv (received_signal, Fc, tau, Beta, Fs)

% This function : Receved acquisition signal in real time
% calculate reference function, and pulse compression
% If the acquision signal dos't exis reference function is used. 

% Specs: tau=pulse width, Beta=frequency sweep range, 
% Fs = sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    persistent count;
    persistent h; 
    persistent x;

    %tau     = 5e-6; Beta=15e6; Fs=40e6;  
    Ts      = 1/Fs;

    % Init stage for the first entry

    %Fc processing
    
    nyquistZone = floor(Fc /(Fs/2));   
    %For Nyquist Zone = 2, 4, 6, ... chirp must be reversed
    if mod(nyquistZone, 2) == 1
        Beta = -Beta;
    end
       
    count   = 1;    
    N       = ceil(tau*Fs);  %Number of sample points of input chirp signal
    t       = 0:Ts:Ts*(N-1);
    c       = 3e8;
    wf_type = 1;   % waveform type 1-Chirp, 2-LFM

    if wf_type==1
        % LFM complex IQ signal    
        mu  = Beta/tau;
        t0  = abs(tau/2);
        x   = exp(1j*pi*mu*((t-t0).^2));    
    else
        % Chirp IF     
        x  = chirp(t,0,tau,Beta);
    end    

    % Multiply reference by Hamming
    x   = x.*hamming(N)';

    % Plot Reference signal
    %figure; plot(t*1e6,real(x)); 
    %title('Transmitted LFM Chirp Pulse'); xlabel('Time (us)');grid

    % Matched Digital FIR filter to LFM Chirp Signal x(t)
    h = ifft(conj(fft(x)));
    
    h = IqModulator2 (h, Fc, 0, Fs);

    % No Input acqusition signal
    if nargin == 0
        received_signal =x;
    end


    % Process Received Signal Through Filter 
    filter_out = conv(h,received_signal);
end

function [filter_out]=r_conv2 (received_signal, h)

% This function : Receved acquisition signal in real time
% calculate reference function, and pulse compression
% If the acquision signal dos't exis reference function is used. 

% Specs: tau=pulse width, Beta=frequency sweep range, 
% Fs = sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process Received Signal Through Filter 
    filter_out = conv(h,received_signal);
end

function [filter_out]=GetMatchedFilter (Fc, tau, Beta, Fs)

% This function : Receved acquisition signal in real time
% calculate reference function, and pulse compression
% If the acquision signal dos't exis reference function is used. 

% Specs: tau=pulse width, Beta=frequency sweep range, 
% Fs = sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    persistent count;
    persistent h; 
    persistent x;

    %tau     = 5e-6; Beta=15e6; Fs=40e6;  
    Ts      = 1/Fs;

    % Init stage for the first entry

    %Fc processing
    
    nyquistZone = floor(Fc /(Fs/2));   
    %For Nyquist Zone = 2, 4, 6, ... chirp must be reversed
    if mod(nyquistZone, 2) == 1
        Beta = -Beta;
    end
       
    count   = 1;    
    N       = ceil(tau*Fs);  %Number of sample points of input chirp signal
    t       = 0:Ts:Ts*(N-1);
    c       = 3e8;
    wf_type = 1;   % waveform type 1-Chirp, 2-LFM

    if wf_type==1
        % LFM complex IQ signal    
        mu  = Beta/tau;
        t0  = abs(tau/2);
        x   = exp(1j*pi*mu*((t-t0).^2));    
    else
        % Chirp IF     
        x  = chirp(t,0,tau,Beta);
    end    

    % Multiply reference by Hamming
    x   = x.*hamming(N)';

    % Plot Reference signal
    %figure; plot(t*1e6,real(x)); 
    %title('Transmitted LFM Chirp Pulse'); xlabel('Time (us)');grid

    % Matched Digital FIR filter to LFM Chirp Signal x(t)
    h = ifft(conj(fft(x)));
    
    h = IqModulator2 (h, Fc, 0, Fs);

    % Process Received Signal Through Filter 
    filter_out = h;
end
%% Author: Joan <Joan@ARBITRARYLAPTOP>
%% Created: 2020-01-21

%% It returns an IQ modulated carrier from the complex baseband signal

function retval = IqModulator2 (iqWfm, cFreq, cPhase, sRate)
  
  % The carrier complex waveform is calculated
  carrier = 0:(length(iqWfm) - 1);
  carrier = carrier / sRate;
  %
  %carrier = cos(2*pi*cFreq * carrier) + 1i * sin(2*pi*cFreq * carrier);
  carrier = exp(1i * 2 * pi * cFreq * carrier +  cPhase);
  % IQ modulation is performed by complex multiplication of the complex baseband
  % and the complex carrier
  retval = iqWfm .* carrier;
  % The actual modulated waveform is just the real part of the product
  retval = real(retval);

end

function [mWfmOut] = ApplyFreqOffset2(fOffset, sampleRate, myWfm)    
    
    wfmL = length(myWfm);
    angleArray = 0:(wfmL - 1);
    angleArray = 2 * pi * fOffset * angleArray;
    angleArray = angleArray / sampleRate;
    
    angleArray = exp(1i * angleArray);
    
    mWfmOut = myWfm .* angleArray;
      
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



