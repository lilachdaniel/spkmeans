#!/bin/bash
#Script to compile and execute a C program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c spkmeans_utils.c -lm -o spkmeans