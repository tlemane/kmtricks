#include <criterion/criterion.h>
#include "code.hpp"

using namespace km;
typedef uint64_t kt;

Test(code, code_encode_uchar)
{
Code<kt> code;
char a = 'a';
char A = 'A';
char t = 't';
char T = 'T';
char c = 'c';
char C = 'C';
char g = 'g';
char G = 'G';
cr_assert(code.encode(a) == 0);
cr_assert(code.encode(A) == 0);
cr_assert(code.encode(t) == 2);
cr_assert(code.encode(T) == 2);
cr_assert(code.encode(c) == 1);
cr_assert(code.encode(C) == 1);
cr_assert(code.encode(g) == 3);
cr_assert(code.encode(G) == 3);
}

Test(code, code_encode_string)
{
Code<kt> code;
kt val = code.encode("ACGTACGT", 8);
cr_assert(val == 0x1E1E);
}

Test(code, code_decode_uchar)
{
Code<kt> code;
cr_assert(code.decode(0) == "AAAA");
cr_assert(code.decode(1) == "AAAC");
cr_assert(code.decode(2) == "AAAT");
cr_assert(code.decode(3) == "AAAG");
}

Test(code, code_decode_value)
{
Code<kt> code;
cr_assert(code.decode(0x1E1E, 8) == "ACGTACGT");
}

Test(code, code_set_custom_encoding)
{
Code<kt> code;
kt val = code.encode("ACGTACGT", 8);
cr_assert(val == 0x1E1E);

uchar newmap[4] = {'T', 'A', 'C', 'G'};
code.set_encoding(newmap);
val = code.encode("ACGTACGT", 8);
cr_assert(val == 0x6C6C);
}