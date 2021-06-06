% =========================================================================
% Copyright (C) 2016-2021 Tabor-Electronics Ltd <http://www.taborelec.com/>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>. 
% =========================================================================
% Author: Nadav Manos, Fractions by Joan Mercade
% Date: May 17, 2021
% Version: 2.0.1


classdef TEProteusInst < handle
    % TEProteusInst: NI-VISA based connection to Proteus Instrument.
    
   
    properties
        ParanoiaLevel = 1; % Paranoia level (0:low, 1:normal, 2:high)        
    end
    
    properties (SetAccess=private)
        ConnStr = ''; % The Connection-String        
        ViSessn = 0;  % VISA Session
    end
    
    properties (Constant=true)
        VISA_IN_BUFF_SIZE = 819200;   % VISA Input-Buffer Size (bytes)
        VISA_IN_BUFF_SIZE_LONG = 8192000;   % VISA Input-Buffer Size for Long Transfers (bytes)
        VISA_OUT_BUFF_SIZE = 819200;  % VISA Output-Buffer Size (bytes)
        VISA_OUT_BUFF_SIZE_LONG = 8192000;  % VISA Output-Buffer Size for Long Transfers (bytes)
        VISA_TIMEOUT_SECONDS = 10;  % VISA Timeout (seconds)
        BINARY_CHUNK_SIZE = 409600;   % Binary-Data Write Chunk Size (samples)
        WAIT_PAUSE_SEC = 0.02;      % Waiting pause (seconds)
    end
    
    methods % public
        
        function obj = TEProteusInst(connStr, paranoiaLevel)
            % TEProteusInst - Handle Class Constructor
            %
            % Synopsis
            %   obj = TEProteusInst(connStr, [verifyLevel])
            %
            % Description
            %   This is the constructor of the VisaConn (handle) class.
            %
            % Inputs ([]s are optional)
            %   (string) connStr      connection string: either a full  
            %                         VISA resource name, or an IP-Address.
            %   (int) [paranoiaLevel = 1] paranoia level [0,1 or 2].
            % 
            % Outputs
            %   (class) obj      VisaConn class (handle) object.
            %
            
            assert(nargin == 1 || nargin == 2);
            
            ipv4 = '^(?:[0-9]{1,3}\.){3}[0-9]{1,3}$';
            if 1 == regexp(connStr, ipv4)
                connStr = sprintf('TCPIP0::%s::5025::SOCKET', connStr);
            end
            
            if nargin == 2
                %verifyLevel = varargin(1);
                if paranoiaLevel < 1
                    obj.ParanoiaLevel = 0;
                elseif paranoiaLevel > 2
                    obj.ParanoiaLevel = 2;
                else
                    obj.ParanoiaLevel = fix(paranoiaLevel);
                end
            else
                obj.ParanoiaLevel = 1;
            end
            
            obj.ConnStr = connStr;
            % Select the right one for the active VISA Library
            obj.ViSessn = visa('NI', connStr);
            %obj.ViSessn = visa('keysight', connStr);
            %obj.ViSessn = visa('tek', connStr);
            
            set(obj.ViSessn, 'OutputBufferSize', obj.VISA_OUT_BUFF_SIZE);
            set(obj.ViSessn, 'InputBufferSize', obj.VISA_IN_BUFF_SIZE);
            obj.ViSessn.Timeout = obj.VISA_TIMEOUT_SECONDS;
            %obj.ViSessn.Terminator = newline;
           
        end
        
        function delete(obj)
            % delete - Handle Class Destructor
            %
            % Synopsis
            %   obj.delete()
            %
            % Description
            %   This is the destructor of the VisaConn (handle) class.
            %   (to be called on a VisaConn class object).
            %
                      
            obj.Disconnect();
            delete(obj.ViSessn);
            obj.ViSessn = 0;
        end
        
        function ok = Connect(obj)
            % Connect - open connection to remote instrument.
            %
            % Synopsis
            %    ok = obj.Connect()
            %
            % Description
            %    Open connection to the remote instrument
            %
            % Outputs
            %    (boolean) ok   true if succeeded; otherwise false.
            %
                        
            ok = false;
            try
                if strcmp(obj.ViSessn.Status, 'open')
                    ok = true;
                else
                    fopen(obj.ViSessn);
                    pause(obj.WAIT_PAUSE_SEC);
                    ok = strcmp(obj.ViSessn.Status, 'open');                    
                end                
            catch ex
                msgString = getReport(ex);
                warning('fopen failed:\n%s',msgString);
            end
        end
		
		function Disconnect(obj)
            % Disconnect - close connection to remote instrument.
            %
            % Synopsis
            %   obj.Disconnect()
            %
            % Description
            %    Close connection to remote-instrument (if open).
            
            if strcmp(obj.ViSessn.Status, 'open')
                stopasync(obj.ViSessn);
                flushinput(obj.ViSessn);
                flushoutput(obj.ViSessn);
                fclose(obj.ViSessn);
            end
        end
        
        function [errNb, errDesc] = QuerySysErr(obj, bSendCls)
            % QuerySysErr - Query System Error from the remote instrument
            %
            % Synopsis
            %   [errNb, [errDesc]] = obj.QuerySysErr([bSendCls])
            %
            % Description
            %   Query the last system error from the remote instrument,
            %   And optionally clear the instrument's errors list.
            %
            % Inputs ([]s are optional)
            %   (bool) [bSendCls = false]  
            %           should clear the instrument's errors-list?
            %
            % Outputs ([]s are optional)
            %   (scalar) errNb     error number (zero for no error).
            %   (string) [errDesc] error description.
            
            if ~exist('bSendCls', 'var')
                bSendCls = false;
            end
            
            obj.waitTransferComplete();
            [answer, count, errmsg] = query(obj.ViSessn, 'SYST:ERR?');
            obj.waitTransferComplete();
                        
            if ~isempty(errmsg)
                error('getError() failed: %s', errmsg);
            end
            
            sep = find(answer == ',');
            if (isempty(sep) || count <= 0 || answer(count) ~= newline)
                warning('querySysErr() received invalid answer: "%s"', answer);
                flushinput(obj.ViSessn);
            end
            
            if ~isempty(sep) && isempty(errmsg)
                errNb = str2double(answer(1:sep(1) - 1));
                errmsg = answer(sep(1):end);
                if 0 ~= errNb && nargin > 1 && bSendCls
                    query(obj.ViSessn, '*CLS; *OPC?');
                end
            else
                errNb =  -1;
                if isempty(errmsg)
                    errmsg = answer;
                end               
            end
            
            if nargout > 1
                errDesc = errmsg;
            end
        end       
        
        
        function SendCmd(obj, cmdFmt, varargin)
            % SendCmd - Send SCPI Command to instrument
            %
            % Synopsis
            %   obj.SendCmd(cmdFmt, ...)
            %
            % Description
            %   Send SCPI Command to the remote instrument.
            %
            % Inputs ([]s are optional)
            %   (string) cmdFmt      command string-format (a la printf).
            %            varargin    arguments for cmdFmt
            obj.waitTransferComplete();
            
            if nargin > 2
                cmdFmt = sprintf(cmdFmt, varargin{1:end});                
            end
            
            resp = '';
            errMsg = '';
            respLen = 0;
            
            if obj.ParanoiaLevel == 0
                fprintf(obj.ViSessn, cmdFmt);
                obj.waitTransferComplete();
            elseif obj.ParanoiaLevel == 1
                cmdFmt = strcat(cmdFmt, ';*OPC?');
                [resp, respLen, errMsg] = query(obj.ViSessn, cmdFmt);
            elseif obj.ParanoiaLevel >= 2
                cmdFmt = strcat(cmdFmt, ';:SYST:ERR?');
                [resp, respLen, errMsg] = query(obj.ViSessn, cmdFmt);
            end
            
            if (obj.ParanoiaLevel > 0 && ~isempty(errMsg))
                error('query(''%s\'') failed\n %s', cmdFmt, errMsg);
            elseif (obj.ParanoiaLevel >= 2 && respLen > 0)
                resp = deblank(resp);
                sep = find(resp == ',');
                if ~isempty(sep)
                    errNb = str2double(resp(1:sep(1) - 1));
                    if 0 ~= errNb
                        query(obj.ViSessn, '*CLS; *OPC?');
                        warning('System Error #%d after ''%s'' (%s).', ...
                            errNb, cmdFmt, resp);
                    end
                end
            end
        end
        
        function resp = SendQuery(obj, qformat, varargin)
            % SendQuery - Send SCPI Query to instrument
            %
            % Synopsis
            %   resp = obj.SendQuery(qformat, ...)
            %
            % Description
            %   Send SCPI Query to the remote instrument,
            %   And return the instrument's response (string).
            %
            % Inputs ([]s are optional)
            %   (string) qformat     query string-format (a la printf).
            %            varargin    arguments for qformat
            %
            % Outputs ([]s are optional)
            %   (string) resp     the instrument's response.
            
            obj.waitTransferComplete();
            if nargin == 2
                [resp, respLen, errMsg] = query(obj.ViSessn, qformat);
            elseif nargin > 2
                qformat = sprintf(qformat, varargin{1:end});
                [resp, respLen, errMsg] = query(obj.ViSessn, qformat);
            else
                resp = '';
                errMsg = '';
                respLen = 0;
            end
            
            if ~isempty(errMsg)
                error('query(''%s\'') failed\n %s', qformat, errMsg);
            end
            
            if respLen > 0
                % remove trailing blanks
                resp = deblank(resp);
            end
        end
        
        function SendBinaryData(obj, pref, datArray, elemType)            
            % SendBinaryData - Send binary data to instrument
            %
            % Synopsis
            %   obj.SendBinaryData(pref, datArray, elemType)
            %
            % Description
            %   Send array of basic-type elements to the remote instrument
            %   as binary-data with binary-data header and (optional) SCPI
            %   statement prefix (e.g. ":TRAC:DATA").
            %
            % Inputs ([]s are optional)
            %   (string) pref      SCPI statement (e.g. ":TRAC:DATA")
            %                      sent before the binary-data header.
            %   (array)  datArray  array of fixed-size elements.
            %   (string) elemType  element type name (e.g. 'uint8')
            
            obj.waitTransferComplete();
            
            
                        
            if ~exist('pref', 'var')
                pref = '';
            end            
            if ~exist('datArray', 'var')
                datArray = [];
            end            
            if ~exist('elemType', 'var')
                elemType = 'uint8';
                datArray = typecast(datArray, 'uint8');
            end 
            
            numItems = length(datArray);  
            switch elemType
                case { 'int8', 'uint8' 'char' }
                    itemSz = 1;
                case { 'int16', 'uint16' }
                    itemSz = 2;
                case { 'int32', 'uint32', 'single' }
                    itemSz = 4;
                case { 'int64', 'uint64', 'double' }
                    itemSz = 8;
                otherwise
                    error('unsopported element-type ''%s''', elemType);
            end
            
            assert(itemSz >= 1 && itemSz <= obj.BINARY_CHUNK_SIZE);
            
            getChunk = @(offs, len) datArray(offs + 1 : offs + len);
            
            % make binary-data header
            szStr = sprintf('%lu', numItems * itemSz);
            pref = sprintf('*OPC?;%s#%u%s', pref, length(szStr), szStr);
            % send it (without terminating new-line!):            
            fwrite(obj.ViSessn, pref, 'char');
            obj.waitTransferComplete();
            
            % send the binary-data (in chunks):            
            offset = 0;
            chunkLen = fix(obj.BINARY_CHUNK_SIZE / itemSz);
            while offset < numItems
                if offset + chunkLen > numItems
                    chunkLen = numItems - offset;
                end
                dat = getChunk(offset, chunkLen);
                fwrite(obj.ViSessn, dat, elemType);
                obj.waitTransferComplete();                
                offset = offset + chunkLen;
            end
            
            % read back the response to that *OPC? query:
            q = fscanf(obj.ViSessn, '%s');
            %fgets(obj.ViSessn, 2);
            
            if obj.ParanoiaLevel >= 2
                [errNb, errDesc] = obj.QuerySysErr(1);
                if 0 ~= errNb
                    warning('System Error #%d (%s) after sending ''%s ..''.', errNb, errDesc, pref);
                end
            end
        end
        
        function datArray = ReadBinaryData(obj, pref, elemType)            
            % ReadBinaryData - Read binary data from instrument
            %
            % Synopsis
            %   datArray = obj.ReadBinaryData(pref, elemType)
            %
            % Description
            %   Read array of basic-type elements from the instrument.
            %
            % Inputs ([]s are optional)
            %   (string) pref      SCPI statement (e.g. ":TRAC:DATA")
            %                      sent before the binary-data header.
            %   (string) elemType  element type name (e.g. 'uint8')
            %
            % Outputs ([]s are optional)
            %   (array)  datArray  array of fixed-size elements.
            
            obj.waitTransferComplete();
            
            set(obj.ViSessn, 'InputBufferSize', obj.VISA_IN_BUFF_SIZE_LONG);
            
            if ~exist('pref', 'var')
                pref = '';
            end            
            
            switch elemType
                case { 'int8', 'uint8' 'char' }
                    itemSz = 1;
                case { 'int16', 'uint16' }
                    itemSz = 2;
                case { 'int32', 'uint32', 'single' }
                    itemSz = 4;
                case { 'int64', 'uint64', 'double' }
                    itemSz = 8;
                otherwise
                    error('unsopported element-type ''%s''', elemType);
            end
            
            assert(itemSz >= 1 && itemSz <= obj.BINARY_CHUNK_SIZE);            
            
            % Send the prefix (if it is not empty)
            if ~isempty(pref)
                fprintf(obj.ViSessn, pref);
            end
            obj.waitTransferComplete();
            
            % Read binary header
            while true
                ch = fread(obj.ViSessn, 1, 'char');
                if ch == '#'
                    break
                end
            end
            
            % Read the first digit
            ch = fread(obj.ViSessn, 1, 'char');
            assert ('0' < ch && ch <= '9');
            
            ndigits = ch - '0';
            %fprintf('ReadBinaryData: ndigits = %d\n', ndigits);
            
            sizestr = fread(obj.ViSessn, ndigits, 'char');
            numbytes = 0;
            for n = 1:ndigits
                ch = sizestr(n, 1);
                numbytes = numbytes * 10 + (ch - '0');
            end
            
            %fprintf('ReadBinaryData: numbytes = %d\n', numbytes);
            
            datLen = ceil(numbytes / itemSz);
            assert(datLen * itemSz == numbytes);
            datArray = zeros(1, datLen, elemType);
            
            chunkLen = fix(obj.BINARY_CHUNK_SIZE / itemSz);
            
            %fprintf('ReadBinaryData: datLen=%d, chunkLen=%d\n', datLen, chunkLen);
            
            % send the binary-data (in chunks):            
            offset = 0;
            
            while offset < datLen
                if datLen - offset < chunkLen
                    chunkLen = datLen - offset;
                end
                datArray(offset + 1 : offset + chunkLen) = ...
                    fread(obj.ViSessn, chunkLen, elemType);
                %obj.waitTransferComplete();                
                offset = offset + chunkLen;
            end
            
            % read the terminating newline character
            ch = fread(obj.ViSessn, 1, 'char');
            assert(ch == newline);
            
            set(obj.ViSessn, 'InputBufferSize', obj.VISA_IN_BUFF_SIZE);
        end    
        
        function model = identifyModel(obj)
            idnStr = obj.SendQuery('*IDN?');    
            idnStr = split(idnStr, ',');    

            if length(idnStr) > 1
                model = idnStr(2);
            else
                model ='';
            end

            model = char(model);
        end

        function options = getOptions(obj)
            optStr = obj.SendQuery('*OPT?');    
            options = split(optStr, ',');    
        end

        function maxSr = getMaxSamplingRate2(obj, model)

            maxSr = 9.0E+9;

            if contains(model, 'P258')
                maxSr = 2.5E+9;
            elseif contains(model, 'P128')
                maxSr = 1.25E+9;
            end
        end
        
        function maxSr = getMaxSamplingRate(obj)
            maxSr = obj.SendQuery(':FREQ:RAST MAX?');    
            maxSr = str2double(maxSr);
        end

        function minSr = getMinSamplingRate2(obj, model)

            minSr = 1.0E+9;    
        end
        
        function minSr = getMinSamplingRate(obj)
            minSr = obj.SendQuery(':FREQ:RAST MIN?');    
            minSr = str2double(minSr);
        end

        function granularity = getGranularity(obj, model, options)

            flagLowGranularity = false;

            for i = 1:length(options)
                if contains(options(i), 'LWG')
                    flagLowGranularity = true;
                end        
            end
            
            sR = obj.SendQuery(':FREQ:RAST?');    
            sR = str2double(sR);
            % For P9082 and P9484 granularity is 64 for SR > 2.5E9
            granularity = 64;    
            if flagLowGranularity && sR<=2.5E9
                granularity = 32;
            end        

            if contains(model, 'P258')
                granularity = 32;    
                if flagLowGranularity
                    granularity = 16;
                end
            elseif contains(model, 'P128')
                granularity = 32;    
                if flagLowGranularity
                    granularity = 16;
                end
            end
        end

        function numOfChannels = getNumOfChannels(obj, model)

            numOfChannels = 4;

            if contains(model, 'P9082')
                numOfChannels = 2;
            elseif contains(model, 'P9482')
                numOfChannels = 2;
            elseif contains(model, 'P1282')
                numOfChannels = 2;
            elseif contains(model, 'P2582')
                numOfChannels = 2;
           elseif contains(model, 'P2584')
                numOfChannels = 4;
            elseif contains(model, 'P2588')
                numOfChannels = 8;
            elseif contains(model, 'P25812')
                numOfChannels = 12;
            end
        end
        
        function dacRes = getDacResolution2(obj, model)

            dacRes = 16;

            if contains(model, 'P908')
                dacRes = 8;            
            end
        end
        
        function dacRes = getDacResolution(obj)           

            dacRes = obj.SendQuery(':TRAC:FORM?'); 
            
            if contains(dacRes, 'U8')
                dacRes = 8;
            else
                dacRes = 16;
            end
        end
        
        function retval = Quantization (obj, myArray, dacRes)
  
          minLevel = 0;
          maxLevel = 2 ^ dacRes - 1;  
          numOfLevels = maxLevel - minLevel + 1;

          retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
          retval = retval + minLevel;

          retval(retval > maxLevel) = maxLevel;
          retval(retval < minLevel) = minLevel;

        end       
        
        
    end % public methods
    
    methods (Access = private) % private methods
        
        function waitTransferComplete(obj)
            % waitTransferComplete - wait till transfer status is 'idle'
            while ~strcmp(obj.ViSessn.TransferStatus,'idle')
                pause(obj.WAIT_PAUSE_SEC);
            end
        end
    end % private methods
    
end

