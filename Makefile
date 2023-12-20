all: tests examples

tests:
	cd tests && $(CXX) -DAPSFIND -o apsfind-tests tests.cpp ../apsfind.cpp -I.. && echo "Running ApsFind tests" && ./apsfind-tests > apsfind-tests.out.txt
	cd tests && $(CXX) -DBISECTION -o bisection-tests tests.cpp bisection.cpp -I.. && echo "Running bisection tests for comparison" && ./bisection-tests > bisection-tests.out.txt

examples:
	$(CXX) -c apsfind.cpp
	cd examples && $(CC) -o example_c example.c elliptic_integral.c ../apsfind.o -lm && echo "Running C example" && ./example_c
	cd examples && $(CXX) -o example_cpp example.cpp elliptic_integral.c ../apsfind.o -lm && echo "Running C++ example" && ./example_cpp
	cd examples && $(CC) -o example-uni example-uni.c ../apsfind.o -lm && echo "Running C example for univariate" && ./example-uni

clean:
	rm -f tests/apsfind-tests tests/apsfind-tests.out.txt tests/bisection-tests tests/bisection-tests.out.txt apsfind.o examples/example_c examples/example_cpp examples/example-uni

.PHONY: tests examples clean
