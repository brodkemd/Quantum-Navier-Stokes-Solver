g++ -c -Wall -Werror -fpic c++_backend.cpp
g++ -shared -o libbackend.so c++_backend.o
g++ -o run c++_backend.o