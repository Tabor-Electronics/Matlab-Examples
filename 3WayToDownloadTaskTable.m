%====================================================================
%Copyright (C) 2021 Tabor-Electronics Ltd <http://www.taborelec.com/>
%
% This example is for use with Tabor Electronics Proteus series of 
% instruments. It demonstrates how to prepare and download waveforms, 
% markers and task tables. For it to run it is also required to have 
% the TEProteusInst.m file on the same path 
%
% For any questions or assistance please contact Tabor Electronics
% support at <support@taborelec.com>
%===================================================================

close all;
clear variables;
clear global;
clc;


 stringFileName='C:\Users\Proteus\Downloads\taskTable.bin'; % file path in the unit system

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
    
    %% Open Communication to Instrument
   
    %connStr = input('Enter the instrument''s IP address  ', 's');
     connStr = 'TCPIP::192.168.0.101::5025::SOCKET';
    inst = TEProteusInst(connStr, paranoia_level);
    fprintf('\n');
        
    res = inst.Connect();
    assert (res == true);
    
    % Identify instrument using the standard IEEE-488.2 Command
    idnstr = inst.identifyModel();
    fprintf('\nConnected to: %s\n', idnstr);
    
    %%Initialize instrument
    fprintf('Reset instrument ..\n');
    inst.SendCmd('*CLS; *RST');
 
    % Get options using the standard IEEE-488.2 Command
    optstr = inst.getOptions();

    % Get granularity
    granul = inst.getGranularity(idnstr, optstr);
      
    % Get Number of Channels
    num_channels = inst.getNumOfChannels(idnstr);
    
    % Get Model name
    model_name=inst.SendQuery(':SYST:INF:MODel?');

    % Setup model dependant parameters 
    if strncmp(model_name,'P908',4)
        dac_mode = 8;
        bpp = 1;
        wav_offs_factor = 1;
        num_ddr = num_channels;
        markers_per_chan = 4;
    elseif strncmp(model_name,'P948',4)
        dac_mode = 16;
        bpp = 2;
        wav_offs_factor = 1;
        num_ddr = int(num_channels/2);
        markers_per_chan = 2;
    else
        dac_mode = 16;
        bpp = 2;
        wav_offs_factor = 2;
        num_ddr = fix(num_channels/2);
        markers_per_chan = 2;
    end
     
    max_dac = 2^dac_mode - 1;
    half_dac = floor(max_dac / 2.0);
    
    
    % Sampling Rate you should define by your needs and according to the capability of the device
    
    %Set SCLK OR you can get the maximum, sample rate for the target by inst.getMaxSamplingRate();
    samplingRate = 1.e9;  % 1GHz
    inst.SendCmd('FREQ:RAST %d', samplingRate);  
    
    

    %% Build waveforms and markers:
    
    fprintf('Build wave-data and markers-data for 12 segments ..');
    
    % Set Segnment length
    seglen = 2^24; %  points

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
       
    low_dac_level = 0;
    
    tStart = tic;    
    % Build DC level and ramps with the formula y=a*x+b
    for i=1:num_channels   
       a=0;
       b=0;  
       switch i   
            case 1  % DC level 0%
                b=0;
                a=0;
                x = linspace(1, 0,seglen);
            case 2  % DC level 25%
                b=0.25;
                 a=0;
                x = linspace(1, 0,seglen);
            case 3  % DC level 50%
                b=0.5;
                 a=0;
                x = linspace(1, 0,seglen);
            case 4   % DC level 75%
                b=0.75;
                x = linspace(1, 0,seglen);
            case 5    % DC level 100%
                 b=1;
                 a=0;
                 x = linspace(1, 0,seglen);
            case 6                                
                 a=1;
                 x = linspace(1, 0,seglen);
            case 7 
                b=1;
                a=1;
                 x = linspace(0, 0,seglen);
            case 8
               b=0.5;
               a=1;
                x = linspace(0, 1,seglen);
            case 9
                 b=0.5;
                a=1;
                x = linspace(0, 1,seglen);
            case 10 
                a=1;
                 b=0.25;
                 x = linspace(1, 0,seglen);
            case 11 
                a=1;
                 b=0.50; 
                 x = linspace(1, 0,seglen);
            case 12 
                a=1;
                b=0.25; 
                x = linspace(0, 1,seglen);
       end
             
          y=(a*x+b)*max_dac; 
          
          %plot(y); 
          y=round(y);
          y=min(y, max_dac);
          y=max(y, low_dac_level);
          
       if dac_mode == 16
            waves(i,:)=uint16(y);
       else
            waves(i,:)=uint8(y);
       end   
     end

    fprintf('done preparing data\n')
    
    %% Download waveforms
    for j=1:num_ddr

        if dac_mode == 16
        channb=2*j-1;
        else
        channb=j;
        end

        % Select channel
        inst.SendCmd( sprintf(':INST:CHAN %d', channb));        
        
        for jj=1:num_channels

            segnum = jj;
            wav = waves(jj,:);
            mrk = marks(jj,:);
            
            tic
            fprintf('Download wave to segment %d of channel %d\n', segnum, channb);

            %define segment
            inst.SendCmd( sprintf(':TRAC:DEF %d, %d', segnum, seglen) );

            %Select the segment
            inst.SendCmd( sprintf(':TRAC:SEL %d', segnum) );

            %%%set data type uint8 or uint16
            precision = class(wav);
            
            % Download waveform binary-data in chunks
            
            chunk = 16 * 1024 * 1024; % 16M
            offset = 0;
            totlen = length(wav);
            
            while (offset < totlen)
                chunk = min(chunk, totlen - offset);
                pref = sprintf(':TRAC:DATA %d', offset * wav_offs_factor);
                inst.SendBinaryData(pref, wav, precision);
                offset = offset + chunk;
            end 
            
            out=inst.SendQuery(':SYST:ERR?');
            
           %   Play the specified segment at the selected channel:
            inst.SendCmd(sprintf(':SOUR:FUNC:MODE:SEGM %d', segnum) );
            inst.SendCmd(':SOUR:VOLT 1.0');
            toc

        end       
    end    
          
    fprintf('Finish downloading 12 waves of size %d bytes \n', seglen ); 
    tEnd = toc(tStart)
    
    
    %% There are 3 ways to use task table
    % 
    %  (1)  Create Task Tables and download them using SCPI
    %  (2)  Task table download using Binary downloading  
    %  (3)  Download from existing bin file

    prompt = 'Selet the mode you want to use task table \n (1)  Create Task Tables and download them \n (2)  Download from existing bin file \n (3)   Binary downloading \n';    
    taske_table_mode = input(prompt); 
   
 
 if taske_table_mode ==1      
    %Set task table length
    tasklen = 3 ;    
   
    % Create Task Table 1 on channels 1-4 :    
    % Use the internal task-table composer to prepare task-table rows:    
  for a=1:4
    % Allocate 3 task-table rows    
    inst.SendCmd( sprintf(':INST:CHAN %d', a));   
    
    %task 1 raw 1
    inst.SendCmd(sprintf(':TASK:COMP:LENG %d', tasklen));
    
    %Select Task to define
    inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 1));
    
    %Select Task type SING or part of a sequenc
    inst.SendCmd(':TASK:COMP:TYPE SING');
    
     %Set segment to generate at current task
     inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 5));
         
     % No enabling signal 
     %(don't wait for any signal when entering this task)
     inst.SendCmd(':TASK:COMP:ENAB NONE');   
         
    %Define what will be the next task once task is done    
    inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d', 2));
       
    %task 1 raw 2    
    %Select Task to define
    inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 2));
    
    %Select Task type SING or part of a sequenc
    inst.SendCmd(':TASK:COMP:TYPE SING');
    
    %Set segment to generate at current task
    inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 6));
         
    % No enabling signal 
    %(don't wait for any signal when entering this task)
    inst.SendCmd(':TASK:COMP:ENAB NONE');       
     
    %Define what will be the next task once task is done    
    inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',3))

    %task 1 raw 3
    %Select Task to define
    inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 3));
    
    %Select Task type SING or part of a sequenc
    inst.SendCmd(':TASK:COMP:TYPE SING');
    
     %Set segment to generate at current task
     inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 1));
         
     % No enabling signal 
     %(don't wait for any signal when entering this task)
     inst.SendCmd(':TASK:COMP:ENAB NONE');       
     
    %Define what will be the next task once task is done    
    inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',1))
    
    %Select channel
    inst.SendCmd(sprintf(':INST:CHAN %d', a));
    tic
    %Download the task table
    inst.SendCmd(':TASK:COMP:WRIT');
    fprintf('Downloading Task table of channel %d\n',a); 
    toc
    
  end
    
    %Set SCLK
    inst.SendCmd('FREQ:RAST 1e9'); 
   
    % Create Task Table 2 on channels 5-8      
    % Use the internal task-table composer to prepare task-table rows:    
    % Allocate 3 task-table rows
    for a=1:4     
        inst.SendCmd( sprintf(':INST:CHAN %d', a+4)); 

        %task 1 raw 1
        inst.SendCmd(sprintf(':TASK:COMP:LENG %d', tasklen));

        %Select Task to define
        inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 1));

        %Select Task type SING or part of a sequenc
        inst.SendCmd(':TASK:COMP:TYPE SING');

         %Set segment to generate at current task
         inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 1));

         % No enabling signal 
         %(don't wait for any signal when entering this task)
         inst.SendCmd(':TASK:COMP:ENAB NONE');   

        %Define what will be the next task once task is done    
        inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d', 2));

         %task 1 raw 2
        %Select Task to define
        inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 2));

        %Select Task type SING or part of a sequenc
        inst.SendCmd(':TASK:COMP:TYPE SING');

         %Set segment to generate at current task
         inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 8));

         % No enabling signal 
         %(don't wait for any signal when entering this task)
         inst.SendCmd(':TASK:COMP:ENAB NONE');       

        %Define what will be the next task once task is done    
        inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',3))


        %task 1 raw 3
        %Select Task to define
        inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 3));

        %Select Task type SING or part of a sequenc
        inst.SendCmd(':TASK:COMP:TYPE SING');

         %Set segment to generate at current task
         inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 5));

         % No enabling signal 
         %(don't wait for any signal when entering this task)
         inst.SendCmd(':TASK:COMP:ENAB NONE');       

        %Define what will be the next task once task is done    
        inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',1))

        %Select channel
        inst.SendCmd(sprintf(':INST:CHAN %d', a+4));
        tic
        %Download the task table
        inst.SendCmd(':TASK:COMP:WRIT');
         fprintf('Downloading Task table of channel %d\n',a+4); 
         toc
    end 
    
    %Set SCLK
    inst.SendCmd('FREQ:RAST 1e9'); 
    
    % Create Task Table 5 on channels 9-12 :     
    %Set task table length
    tasklen = 3 ;
    
    % Use the internal task-table composer to prepare task-table rows:
    for a=1:4
        % Allocate 3 task-table rows    
        inst.SendCmd( sprintf(':INST:CHAN %d', a+8));   

        %task 1 raw 1
        inst.SendCmd(sprintf(':TASK:COMP:LENG %d', tasklen));

        %Select Task to define
        inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 1));

        %Select Task type SING or part of a sequenc
        inst.SendCmd(':TASK:COMP:TYPE SING');

         %Set segment to generate at current task
         inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 11));

         % No enabling signal 
         %(don't wait for any signal when entering this task)
         inst.SendCmd(':TASK:COMP:ENAB NONE');   

        %Define what will be the next task once task is done    
        inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d', 2));

         %task 1 raw 2

        %Select Task to define
        inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 2));

        %Select Task type SING or part of a sequenc
        inst.SendCmd(':TASK:COMP:TYPE SING');

         %Set segment to generate at current task
         inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 12));

         % No enabling signal 
         %(don't wait for any signal when entering this task)
         inst.SendCmd(':TASK:COMP:ENAB NONE');       

        %Define what will be the next task once task is done    
        inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',3))

        %task 1 raw 3

        %Select Task to define
        inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 3));

        %Select Task type SING or part of a sequenc
        inst.SendCmd(':TASK:COMP:TYPE SING');

         %Set segment to generate at current task
         inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 1));

         % No enabling signal 
         %(don't wait for any signal when entering this task)
         inst.SendCmd(':TASK:COMP:ENAB NONE');       

        %Define what will be the next task once task is done    
        inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',1))

        %Select channel
        inst.SendCmd(sprintf(':INST:CHAN %d', a+8));
        tic
        %Download the task table
        inst.SendCmd(':TASK:COMP:WRIT');
        fprintf('Downloading Task table of channel %d\n',a+8);
        toc
    end
       
    %Set SCLK
    inst.SendCmd('FREQ:RAST 1e9'); 
 end
 
 if  taske_table_mode == 2
    % Simple Sequence with 4 tasks
    totalTasks = 4;
    mySequence = CreateTaskTable(totalTasks);  

    for taskNumber = 1:totalTasks
        % Assigns segment for task in the sequence 1..numOfSegments
        param = mod(taskNumber - 1, num_channels) + 1;
        mySequence(taskNumber) = SetValueInTask(mySequence(taskNumber),...
                            'segNb', param);  
        % Next Task is t
        param = mod(taskNumber, totalTasks) + 1;
        mySequence(taskNumber) = SetValueInTask(mySequence(taskNumber),...
                            'nextTask1', param);    

        param = taskNumber;
        mySequence(taskNumber) = SetValueInTask(mySequence(taskNumber),...
                                'taskLoopCount', param);

    end
    
    
    for channel=1: num_channels
        tic
        %Select Channel
        inst.SendCmd(sprintf(':INST:CHAN %d', channel));
        % DAC Mode set to 'DIRECT" (Default)
        inst.SendCmd(':FUNC:MODE TASK');

        % Convert task table to binary format for download
        binSequence = TaskTableToBin(mySequence);

        fprintf(1, 'SEQUENCE CREATED!\n');

        fprintf(1, 'DOWNLOADING SEQUENCE!\n');

        prefix = ':TASK:DATA';
        inst.SendBinaryData(prefix, binSequence, 'uint8');

        fprintf(1, 'SEQUENCE CREATED!\n');
        toc
        
        % Select segment for generation
        fprintf(1, 'SETTING AWG OUTPUT\n');
    end

    % Output volatge set to MAX
    inst.SendCmd(':SOUR:VOLT MAX');   
    % Activate outpurt and start generation
    inst.SendCmd(':OUTP ON');
    
 end    
 
 if  taske_table_mode == 3
     
    %% Creating a bin file for task table     
    % the reason for that becuase you can not create a file of task table
    % by your self you need to creat a table task down load it and then
    % save it as a file.
             
    %Set task table length
    tasklen = 3 ;    
     
    % Use the internal task-table composer to prepare task-table rows:
    
    % Allocate 3 task-table rows    
    inst.SendCmd( sprintf(':INST:CHAN %d', 1));   

    %task 1 raw 1
    inst.SendCmd(sprintf(':TASK:COMP:LENG %d', tasklen));

    %Select Task to define
    inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 1));

    %Select Task type SING or part of a sequenc
    inst.SendCmd(':TASK:COMP:TYPE SING');

     %Set segment to generate at current task
     inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 5));

     % No enabling signal 
     %(don't wait for any signal when entering this task)
     inst.SendCmd(':TASK:COMP:ENAB NONE');   

    %Define what will be the next task once task is done    
    inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d', 2));

     %task 1 raw 2

    %Select Task to define
    inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 2));

    %Select Task type SING or part of a sequenc
    inst.SendCmd(':TASK:COMP:TYPE SING');

     %Set segment to generate at current task
     inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 6));

     % No enabling signal 
     %(don't wait for any signal when entering this task)
     inst.SendCmd(':TASK:COMP:ENAB NONE');       

    %Define what will be the next task once task is done    
    inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',3))

    %task 1 raw 3

    %Select Task to define
    inst.SendCmd(sprintf(':TASK:COMP:SEL %d', 3));

    %Select Task type SING or part of a sequenc
    inst.SendCmd(':TASK:COMP:TYPE SING');

    %Set segment to generate at current task
     inst.SendCmd(sprintf(':TASK:COMP:SEGM %d', 1));

    % No enabling signal 
    %(don't wait for any signal when entering this task)
    inst.SendCmd(':TASK:COMP:ENAB NONE');       

    %Define what will be the next task once task is done    
    inst.SendCmd(sprintf(':TASK:COMP:NEXT1 %d',1))

    %Select channel
    inst.SendCmd(sprintf(':INST:CHAN %d', 1));
    %Download the task table
    inst.SendCmd(':TASK:COMP:WRIT');
     
    %Set SCLK
    inst.SendCmd('FREQ:RAST 1.25e9'); 
    
    %save the task table to file
    % Select channel
    inst.SendCmd( sprintf(':INST:CHAN %d', 1));  
    precision= class(stringFileName);
    inst.SendBinaryData(':TASK:FILE',  stringFileName, precision);
    fprintf('saving the task table to file from path: %s',stringFileName); 

    inst.SendCmd(':TASK:FILE:STORE ');
         
    %% Open and upload bin file of the task table  to all channels:
    
   
    for channb=1:num_channels
        % Select channel
        inst.SendCmd( sprintf(':INST:CHAN %d', channb));  
      
        tic 
               
        precision= class(stringFileName);

        inst.SendBinaryData(':TASK:FILE',  stringFileName, precision);
        fprintf('Uploading the file from path: %s \n',stringFileName); 

        inst.SendCmd(':TASK:FILE:LOAD ');
        
        toc
        
        fprintf('Downloading Task table of channel %d  \n', channb); 
        
        %Set SCLK
        inst.SendCmd('FREQ:RAST 1.25e9'); 

        %Turn on channel output
        inst.SendCmd(':OUTP ON');

        %Change to task mode
        inst.SendCmd(':FUNC:MODE TASK');

    end 
    
 end  

 %% Change to Task mode and turn on outputs
    for ii=1:num_channels
        channb=ii;

        %Select channel
        inst.SendCmd(sprintf(':INST:CHAN %d', channb));

        %Turn on channel output
        inst.SendCmd(':OUTP ON');
        
        %Change to task mode
        inst.SendCmd(':FUNC:MODE TASK');
        
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
   
    %%  Set all channels wait for Trigger1 of the master module
    
    input('Press Enter to continue to external triger mode ')
        
    % Write the 3 updated task-table rows to each channel
   for ii=1:num_channels

        inst.SendCmd(sprintf(':INST:CHAN %d', ii));     
        inst.SendCmd(':TRIG:SEL TRG1');        
        inst.SendCmd(':TRIG:LEV 0.2');        
        inst.SendCmd(':TRIG:SOUR:ENAB TRG1');         
        inst.SendCmd(':TRIG:STATE ON');
        inst.SendCmd(':TASK:COMP:SEL 1');
        inst.SendCmd(':TASK:COMP:ENAB TRG1');   
        inst.SendCmd(':TASK:COMP:WRIT');
      
   end    
    %Set SCLK
    inst.SendCmd('FREQ:RAST 1e9'); % this command is important to sync all channels
    
    
    %% Close the session
    inst.Disconnect();

    clear inst

catch ME
    close all;
    rethrow(ME)
end

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

function outTaskTable = CreateTaskTable(n)   
    outTaskTable = GetDefaultTask();
    outTaskTable = repmat(outTaskTable, 1, n);
end


function taskEntry = GetTaskEntry(  segN, ...
                                    nxTask1, ... 
                                    nxTask2, ...
                                    taskLoop, ...
                                    seqLoop, ...
                                    nextDelay, ...
                                    dcVal, ...
                                    idlWvf, ...
                                    enabSig, ...
                                    abortSig, ...
                                    condJump, ...
                                    jumpType,...
                                    state, ...
                                    loopTrig, ...
                                    adcTrig) 
                                
    % The segment number
    taskEntry.segNb =               uint32(segN);
    % The Next-Task for Trigger 1 (zero for end)
    taskEntry.nextTask1 =           uint32(nxTask1);
    % The Next-Task for Trigger 2 (zero for end)
    taskEntry.nextTask2 =           uint32(nxTask2);
    % The task loop count (0:2^20-1)
    taskEntry.taskLoopCount =       uint32(taskLoop);
    % The sequence loop count (0:2^20-1).
    taskEntry.seqLoopCount =        uint32(seqLoop);
    % The delay in clocks before executing the next task.
    taskEntry.nextTaskDelay =       uint16(nextDelay);
    % The DAC value of the idle task DC waveform.  
    taskEntry.taskDcVal =           uint16(dcVal);
    % The behavior during idle-time
    % (0: DC, 1: First Point, 2: Current Segment(?))
    taskEntry.taskIdlWvf =          uint8(idlWvf);
    % The enabling signal
    % (0:None, 1:ExternTrig1, 2:ExternTrig2, 3:InternTrig, 
    %  4:CPU, 5:FeedbackTrig, 6:HW-Ctrl(?))
    taskEntry.taskEnableSig =       uint8(enabSig);
    % The aborting signal
    % (0:None, 1:ExternTrig1, 2:ExternTrig2, 3:InternTrig, 
    %  4:CPU, 5:FeedbackTrig, 6:Any)
    taskEntry.taskAbortSig =        uint8(abortSig);
    % How to decide where to jump
    % 0: Next1, 1: By FBTrig-Value, 2: ExtTrig[1/2]->Next[1/2],
    % 3: NextTaskSel(?), 4: Next Scenario)
    taskEntry.taskCondJumpSel =     uint8(condJump);
    % Task abort jump type
    % (0:Eventually, 1:Immediate)
    taskEntry.taskAbortJumpType =   uint8(jumpType);
    % The task state
    % (0:Single,1:First of sequence, 2:Last of sequence, 3:Inside Sequence)
    taskEntry.taskState =           uint8(state);
    % Enable (1) or disable (0) waiting for trigger on looping.
    taskEntry.taskLoopTrigEn =      uint8(loopTrig);
    % If it's non-zero, gen ADC trigger at the end of the current task.
    taskEntry.genAdcTrigger =       uint8(adcTrig);
end

function taskEntry = GetDefaultTask() 

    taskEntry = GetTaskEntry(   1, ...                                
                                1, ...
                                1, ...
                                1, ...
                                1, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0);
end

function taskEntry = SetValueInTask(inpuTask, fieldName, value) 
                                
    taskEntry = inpuTask;
    fieldName = upper(fieldName);
    
    if strcmp(fieldName, upper('segNb'))
        taskEntry.segNb =               uint32(value);
    elseif strcmp(fieldName, upper('nextTask1'))
        taskEntry.nextTask1 =           uint32(value);
    elseif strcmp(fieldName, upper('nextTask2'))
        taskEntry.nextTask2 =           uint32(value);
    elseif strcmp(fieldName, upper('taskLoopCount'))
        taskEntry.taskLoopCount =       uint32(value);
    elseif strcmp(fieldName, upper('seqLoopCount'))
        taskEntry.seqLoopCount =        uint32(value);
    elseif strcmp(fieldName, upper('nextTaskDelay'))
        taskEntry.nextTaskDelay =       uint16(value);
    elseif strcmp(fieldName, upper('taskDcVal'))
        taskEntry.taskDcVal =           uint16(value);
    elseif strcmp(fieldName, upper('taskIdlWvf'))
        taskEntry.taskIdlWvf =          uint8(value);
    elseif strcmp(fieldName, upper('taskEnableSig'))
        taskEntry.taskEnableSig =       uint8(value);
    elseif strcmp(fieldName, upper('taskAbortSig'))
        taskEntry.taskAbortSig =        uint8(value);
    elseif strcmp(fieldName, upper('taskCondJumpSel'))
        taskEntry.taskCondJumpSel =     uint8(value);
    elseif strcmp(fieldName, upper('taskAbortJumpType'))
        taskEntry.taskAbortJumpType =   uint8(value);
    elseif strcmp(fieldName, upper('taskState'))
        taskEntry.taskState =           uint8(value);
    elseif strcmp(fieldName, upper('taskLoopTrigEn'))
        taskEntry.taskLoopTrigEn =      uint8(value);
    elseif strcmp(fieldName, upper('genAdcTrigger'))
        taskEntry.genAdcTrigger =       uint8(value);
    end
end

function taskTableBin = TaskTableToBin(taskEntry) 
                                
    taskTableBin = [];
   
    for i = 1:length(taskEntry)
        val = taskEntry(i).segNb;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).nextTask1;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).nextTask2;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).taskLoopCount;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).seqLoopCount;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).nextTaskDelay;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).taskDcVal;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).taskIdlWvf;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).taskEnableSig;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).taskAbortSig;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).taskCondJumpSel;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).taskAbortJumpType;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).taskState;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];    
        val = taskEntry(i).taskLoopTrigEn;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).genAdcTrigger;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];         
    end
end