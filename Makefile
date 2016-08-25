all: serial parallel

.PHONY: serial parallel clean unit_tests tester
serial:
	@echo "Compiling the serial code..."
	make -C serial

parallel:
	@echo "Compiling the parallel code..."
	make -C parallel

tester:
	@echo "Compiling serial tester..."
	make -C serial serial_tester
	@echo "Compiling parallel tester..."
	make -C parallel parallel_tester

clean:
	rm -rf *.z *.dSYM *.otf *.thumb *tester
	make -C serial clean
	make -C parallel clean

unit_tests:
	@echo "Serial unit testing..."
	make -C unit_tests test
