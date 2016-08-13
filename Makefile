all: serial parallel

.PHONY: serial parallel clean unit_tests tester
serial:
	@echo "Compiling the serial code..."
	make -C serial

parallel:
	@echo "Compiling the parallel code..."
	make -C parallel

tester:
	@echo "Compiling tester..."
	make -C serial tester
	./tester

clean:
	make -C serial clean
	make -C parallel clean

unit_tests:
	make -C unit_tests test
