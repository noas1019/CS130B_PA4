# Makefile
prog4: prog4.cpp
	g++ -std=c++17 prog4.cpp -o prog4

debug: prog4.cpp
	g++ -ggdb3 -O0 -std=c++17 prog4.cpp -o debug

clean:
	-rm prog4
	-rm debug