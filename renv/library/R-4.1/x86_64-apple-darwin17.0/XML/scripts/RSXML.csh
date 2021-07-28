if(`test -n "-L/usr/local/lib -lxml2 -lz -L/usr/local/lib -llzma -lpthread -liconv -lm"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-L/usr/local/lib -lxml2 -lz -L/usr/local/lib -llzma -lpthread -liconv -lm
else
   setenv LD_LIBRARY_PATH -L/usr/local/lib -lxml2 -lz -L/usr/local/lib -llzma -lpthread -liconv -lm
endif

endif
