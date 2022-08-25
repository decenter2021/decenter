
# Still to improve:
# Variable tudat-bundle-app-dir is not necessary for compilation, we may just append it to the other variable of the tudat-bundle

# -------------- Set paths -----------------------------------------------------
# Matlab dirs
matlab-exe := /usr/local/MATLAB/R2021b/bin/matlab
# Tudat installation specific variables
tudat-bundle-dir := /mnt/nfs/home/lpedroso/tudat-bundle
# Tudat app specific variables
tudat-matlab-feedback-dir := ../src-tudat
cpp-src := simulationConstellation.cpp
matlab-feedback-src := DDELQRConstellationFeedback
matlab-check-control-cycles := 5
# ------------------------------------------------------------------------------

# Current dir for logs
mkfile-dir := $(shell pwd)
date-now := $(shell date +"%Y_%m_%d_%H_%M_%S")
tudat-bundle-app-dir := /mnt/nfs/home/lpedroso/tudat-bundle/tudat/tudat-applications/tudat-app
# tudat-bundle-app-dir := $(tudat-bundle-dir)/mnt/nfs/home/lpedroso/tudat-bundle/tudat/tudat-applications/tudat-app

# High-level targets
all: tudat-app
tudat: tudat-app
#matlab: matlab-check.log

# Run target
run: tudat-app
	rm -f ./output/*.dat
	mkdir -p logs
	xterm -e "$(matlab-exe) -nodisplay -nodesktop -nosplash -batch $(matlab-feedback-src) wait 2>&1 | tee logs/run-matlab-server.log && exit" &
	sleep 10 && ./tudat-app 2>&1 | tee logs/run-tudat-app.log

run-tudat: tudat-app
	mkdir -p logs
	./tudat-app $(mkfile-dir)/output 2>&1 | tee logs/run-tudat-app.log

# Compile Tudat app 
tudat-app: $(cpp-src) $(wildcard *.h)
	rm -rf $(tudat-bundle-app-dir) # Clean build directory
	mkdir -p $(tudat-bundle-app-dir) # Create build directory
	mkdir -p logs
	cp $(cpp-src) *.h $(tudat-bundle-app-dir) # Create source files
	cp $(tudat-matlab-feedback-dir)/*.h $(tudat-bundle-app-dir) # Create source files
	printf "%s%s%s" 'TUDAT_ADD_EXECUTABLE(tudat-app "' $(cpp-src) '" $${Tudat_PROPAGATION_LIBRARIES})' > $(tudat-bundle-app-dir)/CMakeLists.txt # Write CMakeLists
	cd $(tudat-bundle-dir) && bash build.sh 2>&1 | tee $(mkfile-dir)/logs/cmake-tudat-app.log # Compile app
	cp $(tudat-bundle-dir)/build/tudat/bin/tudat-app ./ # Retrieve app

# Check matlab for errors
matlab-check: $(wildcard *.m) $(wildcard *.h)
	$(matlab-exe) -nodisplay -nodesktop -nosplash -r "matlab_check = $(matlab-check-control-cycles); $(matlab-feedback-src); exit;" wait 2>&1 | tee logs/matlab-check.log

# Zip results
output:
	$(matlab-exe) -nodisplay -nodesktop -nosplash -r "path_tudat_app = '$(mkfile-dir)'; run('$(tudat-matlab-feedback-dir)/getTudatOutput.m');exit;" wait 2>&1 | tee logs/get-tudat-output.log
	zip ./output/$(date-now)_output.zip ./output/*.dat ./output/output.mat *.cpp *.h *.m makefile ./logs/*.log
	rm output/*.dat

# Clean
clean: confirmation-clean
	rm -rf $(tudat-bundle-app-dir) # Clean build directory
	rm -f tudat-app # Clean tudat-app
#	rm -f matlab-check.log # Clean matlab-check
	rm -f ./output/*.dat
	rm -f ./output/*.mat
	rm -rf ./logs
	rm -f *.asv


confirmation-clean:
	@echo 'Do you really want to delete all output results, logs, and executables? [y/n]' && read ans && [ $${ans:-n} = y ]

.PHONY: clean confirmation-clean help output

# Help 
help:
	@echo 'all: 		Compile tudat app and perform matlab check'
	@echo 'all: 		Compile tudat app'
	@echo 'tudat:		Compile tudat app'
	@echo 'matlab: 	Perform matlab check'
	@echo 'run: 		Run tudat app and matlab server'
	@echo 'run-tudat: 	Run tudat app'
	@echo 'output: 	Save output results to .mat file and zip raw output files alongside with the source code and logs' 
	@echo 'clean: 		Clean output results, logs, and compiled objects (.zip of past simualtions are not deleted)'


