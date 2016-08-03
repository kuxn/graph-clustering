all: serial parallel

.PHONY: serial parallel clean unit_tests
serial:
	@echo "Compiling the serial code"
	make -C serial

parallel:
	@echo "Compiling the parallel code"
	make -C parallel

clean:
	make -C serial clean
	make -C parallel clean

unit_tests:
	make -C unit_tests test
