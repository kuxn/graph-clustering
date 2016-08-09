all: serial parallel

.PHONY: serial parallel clean bench
serial:
	@echo "Compiling the serial code"
	make -C ../serial
	mv ../main_serial .

parallel:
	@echo "Compiling the parallel code"
	make -C ../parallel
	mv ../main_parallel .

bench:
	make all
	chmod +x ./submit.sh
	./submit.sh

plot:
	chmod +x ./plot.sh
	./plot.sh

clean:
	rm -rf *.z *.dSYM *.otf *.thumb main_serial main_parallel
	rm -rf data times output sbatch* slurm*

