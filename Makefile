all: generate_positions

generate_positions: generate_positions.c
	gcc generate_positions.c -o generate_positions -lm -g
	
clean:
	rm -f generate_positions *.txt
