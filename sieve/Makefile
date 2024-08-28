override LINKFLAGS += -shared -L$(BELLE_TOP_DIR)/x86_64-unknown-linux-gnu/opt/lib/so -L$(BELLE_TOP_DIR)/bin -ltrak  -lcom-coil

override DEBUG =-O3
override COMPILEFLAGS += -I$(BELLE_TOP_DIR)/include -DHEP_SHORT_NAMES  -fpic -pipe
SRC = reco.cc userinfo.cc myutils.cc utility.cc psi.cc # CDCCat.cc CDCStripHit.cc TRAK3.cc TRAKElement.cc TRAKSegment.cc CDCCatHit.cc TRAK4.cc TRAKHitCDC.cc TRAKStrack.cc CDCCatLayer.cc TRAKBox.cc TRAKHitSVD.cc TRAKSuper_CDC.cc CDCClust.\
cc TRAKCStrack.cc TRAKHybrid.cc TRAKTrack.cc CDCClustFinder.cc TRAK.cc TRAKCluster.cc TRAKLayer.cc TRAKWire.cc CDCSector.cc TRAKCylinder.cc  TRAKMaterial.cc CDCStrip.cc        TRAK2.cc        TRAKDSSD.cc      TRAKSection.cc               
# If there are several .cc files to be compiled you must point them above                                                                                                                                                                     

%.o: %.cc %.h
        gcc -c $(COMPILEFLAGS) $(DEBUG) -std=c++11 -MMD $(INC) -o $@ $<

User_reco.so: reco.o userinfo.o myutils.o utility.o psi.o
        gcc $(LINKFLAGS) $(DEBUG)    $^ -o $@
clean:
        -rm *.o *.so *.d


-include $(DEP)

