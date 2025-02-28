/*  PUT - puts command line arguments to stdout, for use in batch files  */

#include <stdio.h>

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        fprintf(stdout, "%s\n", argv[i]);
    }
    return 0;
}
