#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

double parseDouble(char *arg, const char *arg_name)
{
    char *end_ptr;
    double num = strtod(arg, &end_ptr);
    if (errno == ERANGE || end_ptr == arg) {
        fprintf(stderr, "%s must be a valid double-width floating point number, got %s\n", arg_name, arg);
        exit(EXIT_FAILURE);
    }
    return num;
}

unsigned long parseUnsignedLong(char *arg, const char *arg_name)
{
    char *end_ptr;
    unsigned long num = strtoul(arg, &end_ptr, 10);
    if (errno == ERANGE) {
        fprintf(stderr, "value %s for %s is out of range for type unsigned long\n", arg, arg_name);
        exit(EXIT_FAILURE);
    }
    else if (errno == EINVAL || end_ptr == arg) {
        fprintf(stderr, "%s must be a valid positive integer, got %s\n", arg_name, arg);
        exit(EXIT_FAILURE);
    }
    return num;
}
