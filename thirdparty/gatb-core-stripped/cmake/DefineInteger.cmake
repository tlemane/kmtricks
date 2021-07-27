################################################################################
# Got from https://github.com/hoytak/hashreduce/blob/master/cmake/CheckInt128.cmake
################################################################################

include(CheckTypeSize)

MACRO(CHECK_INT128 INT128_NAME VARIABLE DEFINE_NAME)

    if(NOT INT128_FOUND)
        check_type_size("${INT128_NAME}" "int128_t_${DEFINE_NAME}")
        if("HAVE_int128_t_${DEFINE_NAME}")
            if("int128_t_${DEFINE_NAME}" EQUAL 16)
                message("Found: Enabling support for 128 bit integers using ${INT128_NAME}.")
                SET(INT128_FOUND 1)
                SET(${VARIABLE} "${DEFINE_NAME}")
            else()
                message("${INT128_NAME} has size ${int128_t}, can't use.")
            endif()
        endif()
    endif()
endmacro()


################################################################################
# Define INTEGER_KIND and KMER_PRECISION with k as entry
################################################################################

MACRO(DefineInteger K)

    # we initialize the variable
    SET(INT128_FOUND 0)

    Check_Int128 ("__uint128"    128_DEF "USE__uint128")
    Check_Int128 ("__uint128_t"  128_DEF "USE__uint128_t")

    # We may have undefined parameter => we use a default value
    if (NOT k)
        set (k "0")
    endif()

    if (${k} LESS 33)
        set (INTEGER_KIND "1")

    elseif (${k} LESS 65)

        if (INT128_FOUND)
            set (INTEGER_KIND "2")
        else()
            set (INTEGER_KIND "3")
        endif()

    else()
        set (INTEGER_KIND "3")
    endif()


    if (${INTEGER_KIND} EQUAL "3")
        MATH(EXPR KMER_PRECISION "(${k}+31)/32")
    else()
        set (KMER_PRECISION "1")
    endif()

    #message("-- Parameter k=${k} => INTEGER_KIND=${INTEGER_KIND} and KMER_PRECISION=${KMER_PRECISION}")

endmacro()
