mkdir ..\build

g++ -fPIC -c -lm -std=c++17 -Wextra model.cc -o ..\build\model.o
g++ -fPIC -c -lm -std=c++17 -Wextra trace.cc -o ..\build\trace.o
g++ -fPIC -c -lm -std=c++17 -Wextra interptraceclosestpos.cc -o ..\build\interptraceclosestpos.o
