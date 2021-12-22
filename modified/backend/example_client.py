import socket
import struct


# Create a TCP/IP socket
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# Connect the socket to the port where the server is listening
server_address = ('localhost', 9995)
#print('connecting to %s port %s' % server_address)
sock.connect(server_address)

try:
    # Send data
    message = struct.pack('f', 0.2)
    #print('sending "%s"' % message)
    sock.sendall(message)

    data = sock.recv(len(message))
    value = struct.unpack('f', data)
    print('received:', value[0])

finally:
    #print('closing socket')
    sock.close()