%====================================================================
%Copyright (C) 2016 Tabor-Electronics Ltd <http://www.taborelec.com/>
%
%This example is for use with Tabor Electronics Proteus series of 
%instruments. It demonstrates how to prepare and download waveforms, 
%markers and task tables. For it to run it is also required to have 
%the TEProteusInst.m file on the same path 
%
%For any questions or assistance please contact Tabor Electronics
%support at <support@taborelec.com>
%===================================================================

clear all;
close all;
clear variables;
clear global;
clc;

try
% 
%   When sending a SCPI command to instrument, it is highlly recommended to send it as query with *OPC?.
%   For example, instead of sending bare SCPI command :OUTPUT ON, send a compound query :OUTPUT ON; *OPC?
%   In that way the user makes sure the execution of the command has been completed before sending the next command.
%     
%   The file TEProteusInst.m includes a function a called SendCmd that sends the given SCPI string to the instrument.
% 
%   The optional paranoia_level argument can receive the values 0, 1 or 2 where
% 
%   0 means: send scpi_str as a bare SCPI command.
%   1 means: append *OPC? to the given scpi_str, send it as query, and read the response (which is "1").
%   2 means: append :SYST:ERR? to the given, send it as query and read the response (the last system error).
%   If the optional paranoia_level argument is not given then a default paranoia-level (that the user can setup) is used.
% 
%   The initial value of the default paranoia-level is 1.
%   High paranoia-level (2) is good for debugging, because the system-error is checked after each SCPI command.


    paranoia_level = 2; % 0, 1 or 2
    
    %%%%%%%%%Open Communication to Instrument%%%%%%%%%%%%%%%%%%%%
    
    connStr = input('Enter the instrument''s IP address  ', 's');
    %connStr = 'TCPIP::172.16.10.1::5025::SOCKET';
    inst = TEProteusInst(connStr, paranoia_level);
    fprintf('\n');
        
    res = inst.Connect();
    assert (res == true);
    
    %Get instrument ID
    idnstr = inst.SendQuery('*IDN?');
    fprintf('Connected to: %s\n', idnstr);
    
    %%Initialize instrument
    fprintf('Reset instrument ..\n');
    inst.SendCmd('*CLS; *RST');
 
    %get number of channels
    resp=inst.SendQuery(':INST:CHAN? MAX');
    num_channels=str2num(resp)
    
    %Get DAC mode 
    DAC=inst.SendQuery(':SYST:INF:DAC?');

    %set DAC mode and number of DDrs
    if strncmp(DAC,'M0',2)
        dac_mode = 16;
        ddr=num_channels/2;
    elseif strncmp(DAC,'M1',2)
        dac_mode = 8;
        ddr=num_channels
    end

    % Set number of markers per channel
    markers_per_chan = 2;
    
    %Set SCLK 
    inst.SendCmd('FREQ:RAST 2.5e9');

    max_dac = 2^dac_mode - 1;
    half_dac = floor(max_dac / 2.0);


    %%
    % ---------------------------------------------------------------------
    % Build waveforms and markers
    % ---------------------------------------------------------------------
    fprintf('Build wave-data and markers-data for 12 segments ..');
    % Set Segnment length
    seglen = 2^11;
    
    
    num_cycles = [1];
    for n=1:11
        x=2^n;
        num_cycles=[num_cycles,x];
    end

    if dac_mode == 16
        seg_wave_bytes = seglen * 2;
    else
        seg_wave_bytes = seglen;
    end

    seg_mark_bytes = floor(seg_wave_bytes/8);

    if dac_mode == 16
            waves = zeros(12,seglen,'uint16');
        else
            waves = zeros(12,seglen,'uint8');
    end

    marks=zeros(12,seg_mark_bytes, 'uint8');

    num_of_segments = 12;

    for i=1:num_of_segments
        ncycles=num_cycles(i);
        ncycles=1;
        cyclelen=seglen/ncycles;

        if mod(i,3)==0
        %%%%Create Square waveform with varying cycles
            x=0:1:seglen-1;
            y=rem(x, cyclelen);
            y = (y <= cyclelen / 2) * max_dac;
            y=round(y);
            y(y<0)=0;
            y(y>max_dac)=max_dac;

            if dac_mode == 16
                waves(i,:)=uint16(y);
            else
                waves(i,:)=uint8(y);
            end
        end

        if mod(i,3)==1

       %%%%%% build sinusoid wave with varying cycles
            low_dac_level = 0;

            amp = (max_dac - low_dac_level) / 2.0;
            mid = (max_dac + low_dac_level) / 2.0;

            x = linspace(0, 2 * pi, seglen + 1);
            y = sin(ncycles*x(1:seglen)) * amp + mid;

            if dac_mode == 16
                waves(i,:)=uint16(y);
            else
                waves(i,:)=uint8(y);
            end
        end

        if mod(i,3)==2

       %%%%% build triangle wave with varying cycles
           x = linspace(0, 2 * pi, seglen + 1);
           y = sawtooth(ncycles*x(1:seglen)) * amp + mid;

            if dac_mode == 16
                waves(i,:)=uint16(y);
            else
                waves(i,:)=uint8(y);
            end
        end


        if dac_mode == 16
            cycle_bytes = floor(cyclelen/4);
        else
            cycle_bytes = floor(cyclelen/8);
        end
        
        %%%% build markers
        x=0:1:seg_mark_bytes-1;
        y=rem(x, cycle_bytes);
        y = (y <= cycle_bytes / 2) * 255;
        y=round(y);
        y(y<0)=0;
        y(y>max_dac)=max_dac;
        marks(i,:)=uint8(y);
        clear x;
        clear y;

    end

    fprintf('done preparing data\n')
    
    %%%Download waveforms
    for j=1:ddr

        if dac_mode == 16
        channb=2*j-1;
        else
        channb=j;
        end

        % Select channel
        inst.SendCmd( sprintf(':INST:CHAN %d', channb));
        
        
        for jj=1:num_of_segments

            segnum = jj;
            wav = waves(jj,:);
            mrk = marks(jj,:);


            fprintf('Download wave to segment %d of channel %d\n', segnum, channb);

            %define segment
            inst.SendCmd( sprintf(':TRAC:DEF %d, %d', segnum, seglen) );

            %Select the segment
            inst.SendCmd( sprintf(':TRAC:SEL %d', segnum) );

            %%%set data type uint8 or uint16
            precision = class(wav);
            
            %%%% Download the binary data to segment 
            inst.SendBinaryData( ':TRAC:DATA ', wav, precision);
            out=inst.SendQuery(':SYST:ERR?')

    %         if out~=0
    %             fprintf('ERROR: %s after writing binary values\n', out);
    %         end

            fprintf('Download marker to segment %d of channel %d\n', segnum, channb);

            % Download the binary data to marker 
            precision = class(mrk);
            inst.SendBinaryData( ':MARK:DATA', mrk, precision);
            out=inst.SendQuery(':SYST:ERR?')

    %         if out~=0
    %             fprintf('ERROR: %s after writing binary values\n', out);
    %         end

        %   Play the specified segment at the selected channel:
            inst.SendCmd(sprintf(':SOUR:FUNC:MODE:SEGM %d', segnum) );

        end

       
    end

    
    
    
    %%
    %%Create Task Table
    
    %Set task table length
    tasklen = 3;

    for ii=1:num_channels
        channb=ii;
        %Select channel
        inst.SendCmd(sprintf(':INST:CHAN %d', channb));
        
        %Set Task table length
        inst.SendCmd(sprintf(':TASK:COMP:LENG %d', tasklen));

        for jj=1:tasklen
            curr_task = jj;
            loop = jj;
            segnb = jj;
            
            %Select Task to define
            inst.SendCmd(sprintf(':TASK:COMP:SEL %d', curr_task));
            
            %Select Task type SING or part of a sequence
            inst.SendCmd(':TASK:COMP:TYPE SING');
            
            %Set number of loops for current task
            inst.SendCmd(sprintf(':TASK:COMP:LOOP %d', loop));
            
            %Set segment to generate at current task
            inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', segnb));
            
            %Set enable signal to be TRG1. To run the task without trigger
            %comment out this command.
            if curr_task==1
                inst.SendCmd(sprintf(':TASK:COMP:ENAB TRG1'));
            end
            
            %Define what will be the next task once task is done
            if curr_task==tasklen
                inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d', 1));
            else
                inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d', curr_task+1));
            end
        end

        %Download the task table
        inst.SendCmd(':TASK:COMP:WRIT');

        sprintf('Downloading Task table of channel %d',channb)

    end

    %Set trigger level to 0.2V.
    inst.SendCmd(':TRIG:LEV 0.2');
        
    
    %Change to Task mode and turn on outputs
    for ii=1:num_channels

        channb=ii;

        %Select channel
        inst.SendCmd(sprintf(':INST:CHAN %d', channb));

        %Turn on channel output
        inst.SendCmd(':OUTP ON');
        
        %Change to task mode
        inst.SendCmd(':FUNC:MODE TASK');
        
        % Enable TRG1 for current channel. To see an ouptput a valid
        % trigger signal must be connected to TRIG1 IN connector.
        inst.SendCmd(':TRIG:SEL TRG1; :TRIG:STAT ON');
        
        
        %Turn on marker outputs
        if dac_mode==16

            for jj=1:2

                mrk = jj;

                inst.SendCmd(sprintf(':MARK:SEL  %d', mrk));

                inst.SendCmd(':MARK ON');
            end

        elseif dac_mode==8

            for jj=1:4

                mrk = jj;

                inst.SendCmd(sprintf(':MARK:SEL  %d', mrk))

                inst.SendCmd(':MARK ON');
            end
        end
    end
    
         %Set SCLK 
        inst.SendCmd('FREQ:RAST 2.5e9');        
        % Close the session
        inst.Disconnect();

        clear inst

catch ME
    close all;
    rethrow(ME)
end