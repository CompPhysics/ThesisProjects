// Triplet.h
// A sample user-defined data type for illustrating ConfigFile
// Operators << and >> are defined to allow writing to and reading from files
// Richard J. Wagner  24 May 2004

#include <iostream>

struct Triplet
{
	int a, b, c;
	
	Triplet() {}
	Triplet( int u, int v, int w ) : a(u), b(v), c(w) {}
	Triplet( const Triplet& orig ) : a(orig.a), b(orig.b), c(orig.c) {}
	
	Triplet& operator=( const Triplet& orig )
		{ a = orig.a;  b = orig.b;  c = orig.c;  return *this; }
};


std::ostream& operator<<( std::ostream& os, const Triplet& t )
{
	// Save a triplet to os
	os << t.a << " " << t.b << " " << t.c;
	return os;
}


std::istream& operator>>( std::istream& is, Triplet& t )
{
	// Load a triplet from is
	is >> t.a >> t.b >> t.c;
	return is;
}
