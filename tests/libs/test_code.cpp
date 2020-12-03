#include <gtest/gtest.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/code.hpp>

using namespace km;
typedef uint64_t kt;

TEST(code, code_encode_uchar)
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
    EXPECT_EQ(code.encode(a), 0);
    EXPECT_EQ(code.encode(A), 0);
    EXPECT_EQ(code.encode(t), 2);
    EXPECT_EQ(code.encode(T), 2);
    EXPECT_EQ(code.encode(c), 1);
    EXPECT_EQ(code.encode(C), 1);
    EXPECT_EQ(code.encode(g), 3);
    EXPECT_EQ(code.encode(G), 3);
}

TEST(code, code_encode_string)
{
    Code<kt> code;
    kt val = code.encode("ACGTACGT", 8);
    EXPECT_EQ(val , 0x1E1E);
}

TEST(code, code_decode_uchar)
{
    Code<kt> code;
    EXPECT_EQ(code.decode(0), "AAAA");
    EXPECT_EQ(code.decode(1), "AAAC");
    EXPECT_EQ(code.decode(2), "AAAT");
    EXPECT_EQ(code.decode(3), "AAAG");
}

TEST(code, code_decode_value)
{
    Code<kt> code;
    EXPECT_EQ(code.decode(0X1E1E, 8), "ACGTACGT");
}

TEST(code, code_set_custom_encoding)
{
    Code<kt> code;
    kt val = code.encode("ACGTACGT", 8);
    EXPECT_EQ(val, 0x1E1E);

    uchar newmap[4] = {'T', 'A', 'C', 'G'};
    code.set_encoding(newmap);
    val = code.encode("ACGTACGT", 8);
    EXPECT_EQ(val, 0x6C6C);
}