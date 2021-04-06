GITHASH = $(shell git rev-parse --short HEAD)
OPTIONS += -DGITHASH="\"${GITHASH}\""
INCLUDE += $(LOCAL_INCLUDE)

#this file define the var INCLUDE


all: OPTIMISATION = -O2
all: executable

debug: NAME_APPEND =
debug: OPTIMISATION += -g -DDEBUG
debug: executable

profile: NAME_APPEND = _prof
profile: OPTIMISATION += -O2
profile: OPTIONS += -pg
profile: LINK += -pg
profile: executable

executable: $(OBJS)
	$(LINKER) $(LDFLAGS) -o $(EXEC)$(NAME_APPEND) $(OBJS) $(LINK)

-include $(OBJS:.o=.d)

.cpp.o:
	$(COMPILER) $(OPTIMISATION) $(OPTIONS) $(INCLUDE) -c $< > $@
	$(COMPILER) $(DEP) $(OPTIONS) $(INCLUDE) $< > $@.d
	@sed -e 's|.*:|$*.o:|' < $*.o.d > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.o.d | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.o.d
clean:
	rm -f *.o *.d
