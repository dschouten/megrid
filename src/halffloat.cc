/*
  Small floating point class by M Phillips - 2009.
  Based on work by Jeroen van der Zijp, November 2008
  (http://www.fox-toolkit.org/ftp/fasthalffloatconversion.pdf)
  
  1 bit sign
  5 bits exponent
  10 bits significand
  
  This code is provided as is with no warranties or guarantees of
  any kind.
  
*/

#include "halffloat.hh"

#ifdef NO_PRECALCULATION
unsigned int half::mantissatable[2048];
unsigned int half::exponenttable[64];
unsigned short half::offsettable[64];
unsigned short half::basetable[512];
unsigned char half::shifttable[512];

struct halfinit {
  halfinit() {
    half::init();
  }
} dummy;
#endif

namespace std 
{
  half max( half a, half b ) { return ( a < b ? b : a ); }
  half min( half a, half b ) { return ( a > b ? b : a ); }
}

half max( half a, half b ) { return ( a < b ? b : a ); }
half min( half a, half b ) { return ( a > b ? b : a ); }
