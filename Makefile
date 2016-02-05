visibility: visibility.cpp
	clang++ -std=c++11 -O3 -Wall visibility.cpp -o visibility

visibility.dbg: visibility.cpp
	clang++ -std=c++11 -O0 -g -Wall visibility.cpp -o visibility.dbg
