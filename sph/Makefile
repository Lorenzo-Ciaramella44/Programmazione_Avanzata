CC = gcc
CFLAGS = -g -fno-omit-frame-pointer
LM = -lm

HEADERS = sph.h

%.o: %.c Makefile
		   $(CC) $(CFLAGS) -c $< -o $@

info:
	@echo '---------------------------------------------------------------------------------------------------------------------'
	@echo '$$ make testmake 				---> run the program																'
	@echo '$$ make plotvideo				---> create the video from plot														'
	@echo '$$ make clean_o					---> clean up all .o file															'
	@echo '$$ make clean_png				---> clean up all .png file															'
	@echo '$$ make clean_txt				---> clean up all .txt file															'
	@echo '$$ make clean_mp4				---> clean up all .mp4 file															'
	@echo '$$ make info 					---> get make info																	'
	@echo '---------------------------------------------------------------------------------------------------------------------'

testmake: 1sph.o func.o sph.o
	$(CC) $(CFLAGS) -o sph 1sph.o func.o sph.o $(LM)
	@echo '$$ ./sph					---> to run '

plotvideo:
	gnuplot multi_plot.gp
	ffmpeg -framerate 10 -i frame_%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p 1_output.mp4

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

clean_png:
	@if [ ! -f confirmed_clean_png.txt ]; then \
		echo "Do you want to cancel all .png file? (y/n)"; \
		read ans; \
		if [ "$$ans" = "y" ]; then \
			rm -f *.png; \
			touch confirmed_clean_png.txt; \
		else \
			echo "Operation aborted"; \
		fi \
	else \
		rm -f *.png; \
	fi
	rm confirmed_clean_png.txt

clean_mp4:
	echo "Do you want to cancel all .mp4 file? (y/n)"; \
	@rm -i 1_output.mp4
