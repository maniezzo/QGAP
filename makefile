#Makefile working for both Linux and Windows

CXX = g++
CXXFLAGS_WINDOWS = -g -Wall -O -DWINDOWS
CXXFLAGS_LINUX = -g -Wall -O -DLINUX -std=c++0x -fpermissive -m64 -fPIC -fno-strict-aliasing
LDADD_WINDOWS = -lm
LDADD_LINUX = -lcplex -lm -lpthread -ldl

CPX_INCDIR_WINDOWS  = /IBM/ILOG\CPLEX_Studio1263/cplex/include
CPX_LDPATH_WINDOWS  = /IBM/ILOG/CPLEX_Studio1263/cplex/lib/x86_windows_vs2013/stat_mda
CPX_INCDIR_LINUX  = /opt/ibm/ILOG/CPLEX_Studio128/cplex/include/ilcplex
CPX_LDPATH_LINUX  = /opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic
EIGENGITMIRROR_INCDIR_LINUX  = /home/stanekr/git/eigen-git-mirror
SPECTRA_INCDIR_LINUX = /home/stanekr/git/spectra/include

OBJ = QGAP.o 

#%.o: %.cpp
#	$(CXX) $(CXXFLAGS) -I$(CPX_INCDIR_WINDOWS) -c $^ -o $@

DataManager_linux: DataManager.h DataManager.cpp
	$(CXX) $(CXXFLAGS_LINUX) -I$(CPX_INCDIR_LINUX) -I$(EIGENGITMIRROR_INCDIR_LINUX) -I$(SPECTRA_INCDIR_LINUX) -c DataManager.cpp -o DataManager.o

HeuMagnify_linux: HeuMagnify.h HeuMagnify.cpp
	$(CXX) $(CXXFLAGS_LINUX) -I$(CPX_INCDIR_LINUX) -I$(EIGENGITMIRROR_INCDIR_LINUX) -I$(SPECTRA_INCDIR_LINUX) -c HeuMagnify.cpp -o HeuMagnify.o

json_linux: json.h json.cpp
	$(CXX) $(CXXFLAGS_LINUX) -I$(CPX_INCDIR_LINUX) -I$(EIGENGITMIRROR_INCDIR_LINUX) -I$(SPECTRA_INCDIR_LINUX) -c json.cpp -o json.o

QGAP_linux: QGAP.h QGAP.cpp
	$(CXX) $(CXXFLAGS_LINUX) -I$(CPX_INCDIR_LINUX) -I$(EIGENGITMIRROR_INCDIR_LINUX) -I$(SPECTRA_INCDIR_LINUX) -c QGAP.cpp -o QGAP.o

main_linux: main.cpp
	$(CXX) $(CXXFLAGS_LINUX) -I$(CPX_INCDIR_LINUX) -I$(EIGENGITMIRROR_INCDIR_LINUX) -I$(SPECTRA_INCDIR_LINUX) -c main.cpp -o main.o

link_linux: DataManager_linux HeuMagnify_linux json_linux QGAP_linux main_linux
	$(CXX) $(CXXFLAGS_LINUX) -I$(CPX_INCDIR_LINUX) -L$(CPX_LDPATH_LINUX) DataManager.o HeuMagnify.o QGAP.o json.o main.o -o QGAP $(LDADD_LINUX)


DataManager_windows: DataManager.h DataManager.cpp
	$(CXX) $(CXXFLAGS_WINDOWS) -I$(CPX_INCDIR_WINDOWS) -I$(EIGENGITMIRROR_INCDIR_WINDOWS) -I$(SPECTRA_INCDIR_WINDOWS) -c DataManager.cpp -o DataManager.o

HeuMagnify_windows: HeuMagnify.h HeuMagnify.cpp
	$(CXX) $(CXXFLAGS_WINDOWS) -I$(CPX_INCDIR_WINDOWS) -I$(EIGENGITMIRROR_INCDIR_WINDOWS) -I$(SPECTRA_INCDIR_WINDOWS) -c HeuMagnify.cpp -o HeuMagnify.o

json_windows: json.h json.cpp
	$(CXX) $(CXXFLAGS_WINDOWS) -I$(CPX_INCDIR_WINDOWS) -I$(EIGENGITMIRROR_INCDIR_WINDOWS) -I$(SPECTRA_INCDIR_WINDOWS) -c json.cpp -o json.o

QGAP_windows: QGAP.h QGAP.cpp
	$(CXX) $(CXXFLAGS_WINDOWS) -I$(CPX_INCDIR_WINDOWS) -I$(EIGENGITMIRROR_INCDIR_WINDOWS) -I$(SPECTRA_INCDIR_WINDOWS) -c QGAP.cpp -o QGAP.o

main_windows: main.cpp
	$(CXX) $(CXXFLAGS_WINDOWS) -I$(CPX_INCDIR_WINDOWS) -I$(EIGENGITMIRROR_INCDIR_WINDOWS) -I$(SPECTRA_INCDIR_WINDOWS) -c main.cpp -o main.o

main: DataManager_windows HeuMagnify_windows json_windows QGAP_windows main_windows
	$(CXX) $(CXXFLAGS_WINDOWS) -I$(CPX_INCDIR_WINDOWS) -L$(CPX_LDPATH_WINDOWS) DataManager.o HeuMagnify.o QGAP.o json.o main.o -o QGAP.exe $(LDADD_WINDOWS)


clean:
	rm -rf *.o QGAP.exe QGAP 

.PHONY: clean