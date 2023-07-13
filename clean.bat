del build\*.o
del build\coef\*.o
del build\var\*.o
rmdir build\coef
rmdir build\var
rmdir build


rm lib\libjupitermag.dll

cd lib\libspline
call clean.bat
cd ..\..\

cd lib\libcon2020
call clean.bat
cd ..\..\

cd lib\libinternalfield
call clean.bat
cd ..\..\