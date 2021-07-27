// hash.h-- implement compile-time choice of hash function
//
// -DuseJellyHash on the compiler command line will compile the code with
// a hash function from JellyFish.  This is useful for comparing results with
// the implementation of SBT.
//
// Otherwise, the code is compiled with sabuhash as the has function

#ifndef hash_H
#define hash_H

#ifndef useJellyHash
#define useSabuHash
#endif

#ifdef useJellyHash
#define Hash          JellyHash
#define HashCanonical JellyHashCanonical
#define JellyHashSeed 0x006966796C6C656A  // ascii "jellyfi" as little endian
#include "jellyhash.h"
#endif // useJellyHash

#ifdef useSabuHash
#define Hash          SabuHash
#define HashCanonical SabuHashCanonical
#include "sabuhash.h"
#endif // useSabuHash


#endif // hash_H
