# Makefile for calculating power spectrum from mock data
CC 	= g++
CFLAGS	= -fopenmp 
D_UTIL	= ../util/

O_UTIL	= util.o integration.o 
O_ALL	= $(O_UTIL) 

L_FFTW	= -lfftw3 -lfftw3_threads -lpthread
L_RFFTW = -lrfftw -lfftw
L_MPI	= -lfftw3 -lmpi

Recon_DR12_mocks: $(O_ALL) Recon_DR12_mocks.o header_recon.h
	$(CC) $(CFLAGS) -o recon_DR12_mock Recon_DR12_mocks.o $(O_ALL) $FFTW_POST_LINK_OPTS $(L_RFFTW) -lm

# calc_pow_LEG: $(O_ALL) calc_pow_LEG.o header.h
	# $(CC)  $(CFLAGS) -o calc_pow_LEG calc_pow_LEG.o $(O_ALL) $(L_FFTW) -lm

# fit_win: $(O_UTIL) fit_win.o header.h
	# $(CC)  $(CFLAGS) -o fit_win fit_win.o $(O_UTIL) -lm

# calc_cov: $(O_ALL) calc_cov.o header.h
	# $(CC)  $(CFLAGS) -o calc_cov calc_cov.o $(O_ALL) $(L_FFTW) -lm

calc_pow_lowz: $(O_ALL) calc_pow_lowz.o header.h
	$(CC)  $(CFLAGS) -o calc_pow_lowz calc_pow_lowz.o $(O_ALL) $(L_FFTW) -lm

calc_cov_lowz: $(O_ALL) calc_cov_lowz.o header.h
	$(CC)  $(CFLAGS) -o calc_cov_lowz calc_cov_lowz.o $(O_ALL) $(L_FFTW) -lm

clean:
	-rm -f *.o 
