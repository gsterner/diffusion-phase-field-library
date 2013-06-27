::Run this script to build dll
::Run simulation with: python simplot.py

gcc -c diffl.c
gcc -shared -o diffl.dll diffl.o