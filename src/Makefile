LIBNAME=event

OBJS=	 ksurvb.o ksurvg.o inthaz.o survkit.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(SHLIBLD) $(SHLIBLDFLAGS) $(FLIBS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
