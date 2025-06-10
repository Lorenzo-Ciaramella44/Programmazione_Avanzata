CC = gcc
CFLAGS = -g -fno-omit-frame-pointer
LM = -lm

HEADERS = nucleo.h

%.o: %.c Makefile
		   $(CC) $(CFLAGS) -c $< -o $@

info:
	@echo '---------------------------------------------------------------------------------------------------------------------'
	@echo '$$ make testmake 				---> run the program																'
	@echo '$$ make plot					---> plot the data																	'
	@echo '$$ make clean_o					---> clean up all file.o															'
	@echo '$$ make clean_txt				---> clean up all file.txt															'
	@echo '$$ make info 					---> get make info																	'
	@echo '---------------------------------------------------------------------------------------------------------------------'
		   
testmake: 1nucleo.o func.o
	$(CC) $(CFLAGS) -o nucleo 1nucleo.o func.o $(LM)
	@echo '$$ ./nucleo 2000					---> to run '

plot:
	gnuplot hist.gp
      
clean_o:
	@if [ ! -f confirmed_clean.txt ]; then \
		echo "Do you want to cancel all .o file? (y/n)"; \
		read ans; \
		if [ "$$ans" = "y" ]; then \
			rm -f *.txt; \
			touch confirmed_clean_o.txt; \
		else \
			echo "Operation aborted"; \
		fi \
	else \
		rm -f *.o; \
	fi
	rm confirmed_clean_o.txt

clean_txt:
	@if [ ! -f confirmed_clean.txt ]; then \
		echo "Do you want to cancel all .txt file? (y/n)"; \
		read ans; \
		if [ "$$ans" = "y" ]; then \
			rm -f *.txt; \
			touch confirmed_clean.txt; \
		else \
			echo "Operation aborted"; \
		fi \
	else \
		rm -f *.txt; \
	fi
	rm confirmed_clean.txt