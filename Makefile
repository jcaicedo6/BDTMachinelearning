FLAGS1 = -O3 -D DEFECT_NO_EXPLICIT_QUALIFICATION
INCLUDES =  
SPECIALFLAGS=-O -fPIC
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
CFLAGS = $(SPECIALFLAGS) -I.
LFLAGS = -O2
RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) -l pthread

FILE=TMVAApply_ids

$(FILE) : $(FILE).o  
	g++ -O3 -o $(FILE).exe $(FILE).o -g $(RLXX) -L -lg2c -lm -lMinuit -lRooFit  -lRooFitCore -lTMVA
$(FILE).o : $(FILE).C
	g++ $(FLAGS1) -o $(FILE).o -c $(FILE).C -g $(INCLUDES) -I $(RCXX)
