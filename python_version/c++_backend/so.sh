g++-11 -c -Wall -Werror -fpic c++_backend.cpp
g++-11 -shared -o libbackend.so c++_backend.o
g++-11 -o run c++_backend.o