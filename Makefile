default: delaunay.exe


delaunay.exe: src/delaunay.cpp include/mesh_t.hpp
	g++ -O3 $< -I./ -I./include -std=c++11 -DNDEBUG -o $@
