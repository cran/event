LIBNAME=event
LD=ld

OBJS=	 ksurvb.o ksurvg.o inthaz.o survkit.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(LD) $(SHLIBLDFLAGS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
