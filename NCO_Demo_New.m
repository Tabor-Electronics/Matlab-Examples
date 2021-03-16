
    clc;
    %parameters
    % AWG Parameters
    connStr = '192.168.1.48'; % your IP here
    paranoia_level = 0; % 0, 1 or 2
    samplingRate = 9E9;  % Up to 9GE9 when 8-bit mode is active
    volt1 = 0.5;

    % NCO settings
    cfr = [1000e6 1000e6 1000e6 1000e6];
    global phase; 
    global refPhase;

    phase = [0.0 0.0 0.0 0.0];

    refPhase = [0.0 0.0 0.0 0.0];

    %phase2 = phase2 + 180;
    apply6db = true;


    %% Load TEPAdmin.dll which is a .Net Assembly
    % The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 
    fprintf('Connecting to Proteus\n');
    asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

    import TaborElec.Proteus.CLI.*
    import TaborElec.Proteus.CLI.Admin.*

    %% Create Administrator
    inst = TEProteusInst(connStr, paranoia_level);
    fprintf('\n');

    res = inst.Connect();
    assert (res == true);

    idnstr = inst.SendQuery('*IDN?');
    fprintf('Connected to: %s\n', idnstr);

    fprintf('Reset instrument ..\n');
    inst.SendCmd('*CLS; *RST');

    % Set sampling rate for AWG as defined in the preamble.
    inst.SendCmd([':SOUR:FREQ:RAST ' num2str(samplingRate)]);

    % ---------------------------------------------------------------------
    % Play selected segment in selected channel and set amplitude
    % ---------------------------------------------------------------------
    for i=1:4
        inst.SendCmd(sprintf(':INST:CHAN %d', i));

        inst.SendCmd(':SOUR:MODE NCO');    

        inst.SendCmd(sprintf(':SOUR:VOLT %d', volt1));


        inst.SendCmd(sprintf(':SOUR:CFR %d', cfr(i)));

        inst.SendCmd(sprintf(':SOUR:PHAS %d', phase(i)));

        if apply6db
            inst.SendCmd(':SOUR:SIXD ON');        
        else
            inst.SendCmd(':SOUR:SIXD OFF');
        end

        inst.SendCmd(':OUTP ON');

    end

    % Set sampling rate for AWG as defined in the preamble.
    inst.SendCmd([':SOUR:FREQ:RAST ' num2str(samplingRate)]);


    % Close the session
    %inst.Disconnect();



    fprintf('\nNCO DEMO COMPLETED\n');

    fig = uifigure('Position',[100 100 500 350]);   
    
    title = uilabel(fig, ...
                    'Text', 'TABOR PROTEUS NCO DEMO',...
                    'FontSize', 24,...
                    'FontWeight', 'bold',...
                    'Position', [80 300 400 60]...
                    );

    %lbl1 = uilabel(fig, ...
    %                'Text', 'CH1',...
    %                'FontWeight', 'bold',...
    %                'Position', [65 275 100 60]...
    %                );

    lbl2 = uilabel(fig, ...
                    'Text', 'CH2',...
                    'FontWeight', 'bold',...
                    'Position', [65 200 100 60]...
                    );

    lbl3 = uilabel(fig, ...
                    'Text', 'CH3',...
                    'FontWeight', 'bold',...
                    'Position', [65 125 100 60]...
                    );

    lbl4 = uilabel(fig, ...
                    'Text', 'CH4',...
                    'FontWeight', 'bold',...
                    'Position', [65 50 100 60]...
                    );

    %sld1 = uislider(fig,...
    %       'Position',[100 300 300 3],...
    %       'Limits', [0 360],...
    %       'ValueChangingFcn',@(sld1,event)setNcoPhase(inst, 1, event));

    sld2 = uislider(fig,...
            'Position',[100 225 300 3],...
            'Limits', [0 360],...
            'ValueChangingFcn',@(sld2,event)setNcoPhase(inst, 2, event));

    sld3 = uislider(fig,...
            'Position',[100 150 300 3],...
            'Limits', [0 360],...
            'ValueChangingFcn',@(sld3,event)setNcoPhase(inst, 3, event));

    sld4 = uislider(fig,...
            'Position',[100 75 300 3],...
            'Limits', [0 360],...
            'ValueChangingFcn',@(sld4,event)setNcoPhase(inst, 4, event)); 
        
    % Create a push button
    btn = uibutton(fig,'push',...
                   'Position',[100, 250, 100, 22],...
                   'Text', 'Set Ref',...
                   'ButtonPushedFcn', @(btn,event) GetNewRef(sld2, sld3, sld4));
               
    btnw = uibutton(fig,'push',...
                   'Position',[300, 250, 100, 22],...
                   'Text', 'Quit',...
                   'ButtonPushedFcn', @(btn,event) QuitProgram(inst));


function GetNewRef(sld2, sld3, sld4)

    global phase; 
    global refPhase;

    refPhase = phase + refPhase;
    refPhase = SetPhaseInRange(refPhase);

    sld2.Value = 0.0;    
    sld3.Value = 0.0;
    sld4.Value = 0.0;

    phase = [0.0 0.0 0.0 0.0];

end

function QuitProgram(instHandle)

    % Close the session
    instHandle.Disconnect();    
    clear inst;
    
    closereq();
end
    
                                                   
function setNcoPhase(instHandle, chNumber, event)

    global phase; 
    global refPhase;
    
    phase(chNumber) = event.Value;

    instHandle.SendCmd(sprintf(':INST:CHAN %d', chNumber));   
    
    instHandle.SendCmd(sprintf(':SOUR:PHAS %d',...
        SetPhaseInRange(phase(chNumber) + refPhase(chNumber))));

end

function outPhase = SetPhaseInRange(inPhase)
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