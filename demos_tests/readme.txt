Test and reference values for test were checked with:

    OS: debian 9.4 64 bit 

    Software versions:
    python3     3.6.3 (all tests ok) / 3.6.4 (all tests ok) / 3.4.2 (moga test gives different result)
    numpy       1.13.3
    deap        1.2.2
    numba       0.36.2
    numexpr     2.6.4
    pyFFTW      0.10.4
    pytest      3.2.5
    
    Ocelot version:
    Ocelot-dev was downloaded 26.05.2018 from https://github.com/ocelot-collab/ocelot



How to run tests

Run all tests
    python3 -m pytest -v
    
Run all tests from specific folder
    python3 -m pytest -v path/

Run all tests from specific file
    python3 -m pytest -v path/filename_test.py

Run single test from specific file
    python3 -m pytest -v path/filename_test.py::test_function
    


How to update reference values
  
Update reference values is possible only for one specific test function per one run
    python3 -m pytest -v path/filename_test.py --update=test_function
