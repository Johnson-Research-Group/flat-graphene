
#package specific stuff
PACKNAME=flatgraphene
TESTDIR=tests

#useful general variables
VENV=test_venv
ACTIVATE=$(VENV)/bin/activate
VENVPIP=$(VENV)/bin/pip
SANDBOX=sbox/sandbox.py

#create/update virtual environment for testing
update-venv:
	@#operations in order:
	@#make virtual environment
	@#activate virtual environment
	@#install package to be tested in editable form (automatically gets dependencies)
	@#deactivate virtual environment
	@#all above directed to shell (necessary)
	@(\
	virtualenv $(VENV); \
	source $(ACTIVATE); \
	$(VENVPIP) install -e .; \
	deactivate; \
	)

#test using virtual environment
test:
	@#operations in order:
	@#activate virtual environment
	@#run all scripts in $(TESTDIR) directory 
	@#deactivate virtual environment
	@#all directed to shell (necessary)
	@(\
	source $(ACTIVATE); \
	for f in $(TESTDIR)/*.py; do python "$$f"; done; \
	deactivate; \
	)

#development sandbox for on the fly testing
sandbox:
	@(\
	source $(ACTIVATE); \
	python $(SANDBOX); \
	deactivate; \
	)

#remove virtual environment
clean:
	@\rm -rf test_venv

