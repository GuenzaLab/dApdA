#!/bin/bash

for i in `seq 1 5`
do
	cp -v param_calc.f08 tmp.f08
	sed -i "s#DUMMY#${i}#g" tmp.f08
	f95 tmp.f08; ./a.out
	rm -rfv tmp.f08
done
exit
