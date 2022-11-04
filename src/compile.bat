call compileobj.bat

g++ -fPIC -lm -std=c++17 -Wextra -shared libjupitermag.cc ..\build\*.o ..\lib\libcon2020\build\*.o ..\lib\libspline\build\*.o ..\lib\libinternalfield\build\*.o -o libjupitermag.dll
