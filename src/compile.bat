g++ -lm -fPIC -std=c++17 -Wextra -O3 trimstring.cc generateheader.cc -o generateheader.exe
call generateheader.exe
del /F /Q generateheader.exe


call compileobj.bat

g++ -fPIC -lm -std=c++17 -Wextra -O3 -shared libjupitermag.cc ..\build\*.o ..\lib\libcon2020\build\*.o ..\lib\libspline\build\*.o ..\lib\libinternalfield\build\*.o -o ..\lib\libjupitermag.dll
