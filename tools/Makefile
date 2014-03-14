
# CONFIG REQUIRED: SET THIS VARIABLE TO THE PATH OF THE INDRI SOURCE CODE
INDRISRC=/Users/mpetri/dev/indri-5.6
# CONFIG REQUIRED: SET THIS TO THE SDSL INSTALL PATH
SDSLPREFIX=/Users/mpetri/
#SDSLPREFIX=../build/external/sdsl-lite/

include $(INDRISRC)/MakeDefns
SHARED=
INCPATH=-I$(INDRISRC)/include $(patsubst %, -I$(INDRISRC)/contrib/%/include, $(DEPENDENCIES)) -I $(SDSLPREFIX)/include/
LIBPATH=-L$(INDRISRC)/obj  $(patsubst %, -L$(INDRISRC)/contrib/%/obj, $(DEPENDENCIES)) -L $(SDSLPREFIX)/lib/
LIBS=-lindri $(patsubst %, -l%, $(DEPENDENCIES)) -ldivsufsort -lsdsl

all: indri
	$(CXX) $(CXXFLAGS) -std=c++11 indri_to_surf.cpp -o indri_to_surf $(INCPATH) $(LIBPATH) $(LIBS) $(CPPLDFLAGS)
	$(CXX) $(CXXFLAGS) -std=c++11 indri_stem_krovetz.cpp -o indri_stem_krovetz $(INCPATH) $(LIBPATH) $(LIBS) $(CPPLDFLAGS)

indri:
	cd  $(INDRISRC); ./configure
	make -C $(INDRISRC) contrib

clean:
	rm -f $(APP)
