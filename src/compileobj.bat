mkdir ..\build

g++ -fPIC -c -lm -std=c++17 -Wextra -O3 model.cc -o ..\build\model.o
g++ -fPIC -c -lm -std=c++17 -Wextra -O3 trace.cc -o ..\build\trace.o
g++ -fPIC -c -lm -std=c++17 -Wextra -O3 interptraceclosestpos.cc -o ..\build\interptraceclosestpos.o
g++ -fPIC -c -lm -std=c++17 -Wextra -O3 footprint.cc -o ..\build\footprint.o
g++ -fPIC -c -lm -std=c++17 -Wextra -O3 coordconv.cc -o ..\build\coordconv.o