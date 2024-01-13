program: conj_grad.c generate_equation.c
	gcc -g -Wall -o program conj_grad.c generate_equation.c -lm


debug: conj_grad.c generate_equation.c
	gcc -g -Wall -D DEBUG1 -o programdb conj_grad.c generate_equation.c -lm

clear: 
	rm program && rm programdb