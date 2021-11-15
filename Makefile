
#package specific stuff
PACKNAME=flatgraphene
TESTDIR=tests

#useful general variables
VENV=test_venv
ACTIVATE=$(VENV)/bin/activate
VENVPIP=$(VENV)/bin/pip

update-venv:
	@#operations in order:
	@#make virtual environment
	@#activate virtual environment
	@#install package to be tested in editable form (automatically gets dependencies)
	@#deactivate virtual environment
	@#all above directed to shell (necessary)
	(\
	virtualenv $(VENV); \
	source $(ACTIVATE); \
	$(VENVPIP) install -e .; \
	deactivate; \
	)
	@#copy tests into virtual environment
	mkdir -p $(VENV)/bin/tests
	cp tests/*.py $(VENV)/bin/tests

test:
	@#operations in order:
	@#activate venv
	@#run all scripts in test directory 
	@#all directed to shell (necessary)
	(\
	source $(ACTIVATE); \
	python $(TESTDIR)/test_shift.py; \
	deactivate; \
	)
	#for f in $(TESTDIR)/*.py; do python "$f"; done; \

clean:
	@\rm -rf test_venv

