function [] = InitServer(server, server_option_0, server_option_1)
  % this function sends the required bytes that configure the server in the manner determined in QNS_InputData.m
  val = int8(server_option_0);
  write(server, val);

  val = int8(server_option_1);
  write(server, val);

  disp("Check Server window for configuration")
end