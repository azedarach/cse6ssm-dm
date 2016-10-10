/**
 * @file read_slha.c
 * @brief contains implementation of SLHA file reading function
 */

#include "read_slha.h"

#include "micromegas.h"
#include "micromegas_aux.h"

int read_slha_file(char* filename)
{
   int error_code;

   error_code = slhaRead(filename, 0);

   return error_code;
}
