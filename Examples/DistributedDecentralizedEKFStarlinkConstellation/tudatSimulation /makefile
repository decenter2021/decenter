# Matlab dirs
matlab-exe := /usr/local/MATLAB/R2021b/bin/matlab

# Tudat installation specific variables
tudat-bundle-dir := /mnt/nfs/home/lpedroso/tudat-bundle
tudat-bundle-app-dir := /mnt/nfs/home/lpedroso/tudat-bundle/tudat/tudat-applications/tudat-app

# Tudat project source specif variables
cpp-src := main.cpp

# Current dir for logs
mkfile-dir := $(shell pwd)

# High-level targets
all: tudat-app matlab-check.log
tudat: tudat-app
matlab: matlab-check.log

# Run target
run: tudat-app matlab-check.log
	mkdir -p logs
	xterm -e "$(matlab-exe) -nodisplay -nodesktop -batch tudatThrustFeedbackMatlabServer wait 2>&1 | tee logs/run-matlab-server.log && exit" &
	sleep 10 && ./tudat-app 2>&1 | tee logs/run-tudat-app.log

run-tudat: tudat-app
	mkdir -p logs
	./tudat-app $(mkfile-dir)/output 2>&1 | tee logs/run-tudat-app.log

# Compile Tudat app 
tudat-app: $(cpp-src) $(wildcard *.h)
	rm -rf $(tudat-bundle-app-dir) # Clean build directory
	mkdir $(tudat-bundle-app-dir) # Create build directory
	mkdir -p logs
	cp $(cpp-src) *.h $(tudat-bundle-app-dir) # Create source files
	printf "%s%s%s" 'TUDAT_ADD_EXECUTABLE(tudat-app "' $(cpp-src) '" $${Tudat_PROPAGATION_LIBRARIES})' > $(tudat-bundle-app-dir)/CMakeLists.txt # Write CMakeLists
	cd $(tudat-bundle-dir) && bash build.sh 2>&1 | tee $(mkfile-dir)/logs/cmake-tudat-app.log # Compile app
	cp $(tudat-bundle-dir)/build/tudat/bin/tudat-app ./ # Retrieve app

# Check matlab for errors
matlab-check.log: $(wildcard *.m) $(wildcard *.h)
	touch matlab-check.log

# Zip results
output:
	$(matlab-exe) -nodisplay -nodesktop -batch getTudatOutput wait 2>&1 | tee logs/get-tudat-output.log
	zip ./output/output.zip ./output/*.dat *.cpp *.h *.m makefile ./logs/*.log matlab-check.log
	rm output/*.dat

# Clean
clean: confirmation-clean
	rm -rf $(tudat-bundle-app-dir) # Clean build directory
	rm -f tudat-app # Clean tudat-app
	rm -f matlab-check.log # Clean matlab-check
	rm -rf ./output
	rm -rf ./logs
	rm -f *.asv

confirmation-clean:
	@echo 'Do you really want to delete all output results, logs, and executables? [y/n]' && read ans && [ $${ans:-n} = y ]

.PHONY: clean confirmation-clean help output

# Help 
help:
	@echo 'all: 		Compile tudat app and perform matlab check'
	@echo 'tudat:		Compile tudat app'
	@echo 'matlab: 	Perform matlab check'
	@echo 'run: 		Run tudat app and matlab server'
	@echo 'run-tudat: 	Run tudat app'
	@echo 'output: 	Save output results to .mat file and zip raw output files alongside with the source code and logs' 
	@echo 'clean: 		Clean output results, logs, and compiled objects'


