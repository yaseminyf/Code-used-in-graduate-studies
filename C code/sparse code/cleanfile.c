// The problem: Write a program to remove the C comments /* ... */ and the
// C++ comments // from a file.
// The file is declared in command line.
//
// usage: cleanfile <filename>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

FILE *infile, *outfile; // infile: input file; outfile: output file


int main(int argc,char *argv[])
{
  
