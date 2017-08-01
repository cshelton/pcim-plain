CXX=g++
RM=rm -f
CPPFLAGS=-std=c++11
LDFLAGS= -lboost_serialization -lpthread
LDLIBS=

SRCS=datainfo.cpp pcim.cpp demo.cpp traj.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: demo

demo: $(OBJS)
	$(CXX) $(LDFLAGS) -o demo $(OBJS) $(LDLIBS)

depend: .depend


.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) .depend

include .depend

