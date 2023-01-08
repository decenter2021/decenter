classdef tudatMatlabServer
    properties
        server;
        cleanUpObj;
        outputFormat = 0; % 0: vector 1:matrix
        SERVER_PORT;
        SERVER_ADDR;
        N_SATS;
        TUDAT_APP_REQUEST_BYTES;
        TUDAT_APP_RESPONSE_BYTES;
    end
    
    methods
        function obj = tudatMatlabServer(port,addr,NSats,outputFormat)
            obj.SERVER_PORT = port;
            obj.SERVER_ADDR = addr;
            obj.outputFormat = outputFormat;
            obj.N_SATS = NSats;
            obj.TUDAT_APP_REQUEST_BYTES = 8*(1+7*NSats);
            obj.TUDAT_APP_RESPONSE_BYTES = 3*8*NSats;
            fprintf("@MATLAB server: Starting MATLAB server.\n");
            obj.server = udp (obj.SERVER_ADDR,obj.SERVER_PORT,'LocalPort',obj.SERVER_PORT);
            obj.server.InputBufferSize = obj.TUDAT_APP_REQUEST_BYTES;
            obj.server.OutputBufferSize = obj.TUDAT_APP_RESPONSE_BYTES;
            obj.server.ByteOrder = 'littleEndian';
            obj.server.Timeout = 100;
            obj.server.Terminator = '';
            obj.server.DatagramTerminateMode = 'on';
             

            fopen(obj.server);
            fprintf("@MATLAB server: Hosted at %s:%d.\n",obj.SERVER_ADDR,obj.SERVER_PORT);
        end
        
        function delete(obj)
            fclose(obj.server);
            return;
        end
        
        
        function waitForClient(obj)
            while (1)
                fprintf('@MATLAB server: Waiting for Tudat client.\n');
                bytesReceived = 0;
                while bytesReceived < obj.TUDAT_APP_REQUEST_BYTES/8
                    [~,count,~,~,datagramport] =  fread(obj.server,[obj.TUDAT_APP_REQUEST_BYTES 1],'double');
                    bytesReceived = bytesReceived + count;
                end
                %Set client address to respond
                obj.server.RemotePort = datagramport;
                if obj.server.ValuesReceived > 0
                    fprintf('@MATLAB server: Client connected.\n');
                    msg = rand(obj.N_SATS*3,1);
                    %fprintf(obj.server,'%f',msg);
                    fwrite(obj.server,msg,'double');
                    %write(obj.server,msg,obj.TUDAT_APP_RESPONSE_BYTES);
                    break;
                end
            end
        end
        
        function [t,x,state] = getRequest(obj)
            bytesReceived = 0;
            data = zeros(obj.TUDAT_APP_REQUEST_BYTES/8,1);
            while bytesReceived < obj.TUDAT_APP_REQUEST_BYTES/8
                [data_aux,count] = fread(obj.server,[obj.TUDAT_APP_REQUEST_BYTES 1],'double');
                data(bytesReceived+1:bytesReceived+count) = data_aux;
                bytesReceived = bytesReceived + count;
            end
            t = data(1);
            state = 0;
            if (t < 0)
                state = 1;
                x = nan;
                t = nan;
                obj.delete();
                return;
            end
            if obj.outputFormat
                x = reshape(data(2:end),[7 obj.N_SATS]);
            else
                x = data(2:end);
            end
  
        end
        
         function sendResponse(obj,u)
            if obj.outputFormat
                u = reshape(u,[3*obj.N_SATS 1]);
            end
            fwrite(obj.server,u,'double')
        end
        
    end
end

