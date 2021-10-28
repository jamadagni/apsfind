all: tests examples

tests:
	cd tests && $(CXX) -DTOMS748 -o toms748-tests tests.cpp ../toms748.cpp -I.. && echo "Running TOMS748 tests" && ./toms748-tests > toms748-tests.out.txt
	cd tests && $(CXX) -DBISECTION -o bisection-tests tests.cpp bisection.cpp -I.. && echo "Running bisection tests for comparison" && ./bisection-tests > bisection-tests.out.txt

examples:
	$(CXX) -c toms748.cpp
	cd examples && $(CC) -o example_c example.c elliptic_integral.c ../toms748.o -lm && echo "Running C example" && ./example_c
	cd examples && $(CXX) -o example_cpp example.cpp elliptic_integral.c ../toms748.o -lm && echo "Running C++ example" && ./example_cpp

clean:
	rm -f tests/toms748-tests tests/toms748-tests.out.txt tests/bisection-tests tests/bisection-tests.out.txt toms748.o examples/example_c examples/example_cpp

.PHONY: tests examples clean
