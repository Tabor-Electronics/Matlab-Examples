% EXAMPLE FOR NCO MODE IN PROTEUS USING VISA
%===================================================
% This example sets the first NCO for all the designated channels
% and the lets the user control the relative phase using slide controls.
% CH1 is used as the reference and phase is always set to 0.
% The "Set Ref" button in the GUI sets then new reference phase after
% aligning all the phases for all the channels.

clc;

% AWG Parameters

% IP Address
connStr = '192.168.1.48'; 
paranoia_level = 0;

% NCO Frequency
ncoFreq = 1E9;
% Boost amplitude for NCO mode
apply6db = true;

% Connection to Proteus
fprintf('Connecting to Proteus\n');

%% Create Handle
inst = TEProteusInst(connStr, paranoia_level);

res = inst.Connect();
assert (res == true);

% Identify instrument using the standard IEEE-488.2 Command
idnstr = inst.identifyModel();
fprintf('\nConnected to: %s\n', idnstr);

% Set sample rate to maximum for target instrument
samplingRate = inst.getMaxSamplingRate();
% If sampling rate lower than 2.5GHz, NCO frequency set to 500MHz
if samplingRate < 2.5E+9
    ncoFreq = 500E+6;
end

% Set number of channels
numOfChan = inst.getNumOfChannels(idnstr);

% NCO settings
% Carrier Frequency for each channel 1..4
cfr = zeros(1, numOfChan);
cfr(1:numOfChan) = ncoFreq;

% Define phase vectors as global variables
% Initial phase and phase reference are set
global phase; 
global refPhase;
phase = zeros(1, numOfChan);
refPhase = zeros(1, numOfChan);

% Reset and initialize instrument
fprintf('Reset instrument ..\n');
inst.SendCmd('*CLS; *RST');

% Set sampling rate for AWG as defined in the preamble.
inst.SendCmd(':SOUR:FREQ:RAST 2.5E9'); % REMOVE
inst.SendCmd(':SOUR:INT X8'); % REMOVE
inst.SendCmd([':SOUR:FREQ:RAST ' num2str(samplingRate)]);

% ---------------------------------------------------------------------
% Play selected segment in selected channel and set amplitude
% ---------------------------------------------------------------------
for i=1:length(cfr)
    % Select channel in a sequence
    inst.SendCmd(sprintf(':INST:CHAN %d', i));
    % Set NCO mode for the current channel
    inst.SendCmd(':SOUR:MODE NCO');  
    % Set ouput voltage to Maximum
    inst.SendCmd(':SOUR:VOLT MAX');
    % Set NCO Frequency for the current channel
    inst.SendCmd(sprintf(':NCO:CFR2 %d', cfr(i))); %CHANGE to CFR1
    % Set NCO Phase for the current channel
    inst.SendCmd(sprintf(':NCO:PHAS2 %d', phase(i))); %CHANGE to PHAS1
    % Apply amplitude boost to NCO
    if apply6db
        inst.SendCmd(':NCO:SIXD2 ON');        
    else
        inst.SendCmd(':NCO:SIXD2 OFF');
    end
    % Switch on output
    inst.SendCmd(':OUTP ON');
end

% Set sampling rate for AWG as defined in the preamble.
inst.SendCmd([':SOUR:FREQ:RAST ' num2str(samplingRate)]);

% GUI creation for the slide controls and buttons

fig = uifigure('Position',[100 100 500 350]);   

title = uilabel(fig, ...
                'Text', 'TABOR PROTEUS NCO DEMO',...
                'FontSize', 24,...
                'FontWeight', 'bold',...
                'Position', [80 300 400 60]...
                );

if (numOfChan > 1)
    lbl2 = uilabel( fig, ...
                    'Text', 'CH2',...
                    'FontWeight', 'bold',...
                    'Position', [65 200 100 60]);
end

if (numOfChan > 2)
    lbl3 = uilabel( fig, ...
                    'Text', 'CH3',...
                    'FontWeight', 'bold',...
                    'Position', [65 125 100 60]);
end

if (numOfChan > 3)
    lbl4 = uilabel( fig, ...
                    'Text', 'CH4',...
                    'FontWeight', 'bold',...
                    'Position', [65 50 100 60]);
end

if (numOfChan > 1)
    sld2 = uislider(fig,...
            'Position',[100 225 300 3],...
            'Limits', [0 360],...
            'ValueChangingFcn',@(sld2,event)setNcoPhase(inst, 2, event));
end

if (numOfChan > 2)
    sld3 = uislider(fig,...
            'Position',[100 150 300 3],...
            'Limits', [0 360],...
            'ValueChangingFcn',@(sld3,event)setNcoPhase(inst, 3, event));
end

if (numOfChan > 3)
    sld4 = uislider(fig,...
            'Position',[100 75 300 3],...
            'Limits', [0 360],...
            'ValueChangingFcn',@(sld4,event)setNcoPhase(inst, 4, event)); 
end

% Create push buttons
btn = uibutton(fig,'push',...
               'Position',[100, 250, 100, 22],...
               'Text', 'Set Ref',...
               'ButtonPushedFcn', @(btn,event) GetNewRef(sld2, sld3, sld4));

btnw = uibutton(fig,'push',...
               'Position',[300, 250, 100, 22],...
               'Text', 'Quit',...
               'ButtonPushedFcn', @(btn,event) QuitProgram(inst));
           
% The GUI is launched and it remains active until the 'quit' button is hit.
% This application is based in the call-bacl mechanism


function GetNewRef(sld2, sld3, sld4)
    % Sets the new reference when all channels are aligned including any
    % differential delay caused by skew, cabling, scope, etc.
    % After hitting the button, the current phase algnment will be used as
    % the 0 degrees reference.

    global phase; 
    global refPhase;

    refPhase = phase + refPhase;
    refPhase = SetPhaseInRange(refPhase);

    if exist('sld2', 'var')
        sld2.Value = 0.0;    
    end
    
    if exist('sld3', 'var')
        sld3.Value = 0.0;
    end
    
    if exist('sld4', 'var')
        sld4.Value = 0.0;
    end
    % New phase control values
    phase = zeros(1, length(phase));

end

function setNcoPhase(instHandle, chNumber, event)
    % The phase of the NCO is updated as specified by the slide value. This
    % callback function is triggered by any change in the slide controls.
    % Phase is only changed for the associated channel to the slide being
    % used. Phase is changed without any disruption in the output of any
    % channel as NCO phase is changed 'on the fly'.

    global phase; 
    global refPhase;
    
    phase(chNumber) = event.Value;

    instHandle.SendCmd(sprintf(':INST:CHAN %d', chNumber));   
    
    instHandle.SendCmd(sprintf(':NCO:PHAS2 %d',...
        SetPhaseInRange(phase(chNumber) + refPhase(chNumber))));

end

function outPhase = SetPhaseInRange(inPhase)

    % All phases are limited to the 0...360 sexagesimal degrees range
    
    outPhase = inPhase; 
    for i=1:length(inPhase)
        while outPhase(i) < 0
            outPhase(i) = outPhase(i) + 360.0;
        end
        
        while outPhase(i) > 360
            outPhase(i) = outPhase(i) - 360.0;
        end
    end
end

function QuitProgram(instHandle)

    % Close the session
    instHandle.Disconnect();    
    clear inst;    
    closereq();
    fprintf('\nNCO DEMO COMPLETED\n');
end
    
                                                   
