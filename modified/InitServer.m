function [] = InitServer(server, server_option_0, server_option_1)
    % this function sends the required bytes that configure the server in the manner determined in QNS_InputData.m
    write(server, int8(server_option_0), int8(server_option_1));

    disp("Check Server window for configuration")
    pause
end