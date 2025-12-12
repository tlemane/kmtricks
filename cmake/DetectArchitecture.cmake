# DetectArchitecture.cmake
# Detects CPU architecture and sets appropriate compile definitions and flags

# Detect architecture based on CMAKE_SYSTEM_PROCESSOR
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)")
    set(ARCH_ARM64 ON CACHE BOOL "ARM64 architecture detected")
    set(ARCH_X86_64 OFF CACHE BOOL "x86_64 architecture not detected")

    # ARM64 guarantees NEON support (ARMv8-A specification)
    add_compile_definitions(ARCH_ARM64=1)

    # Portable ARM64 flags
    set(ARCH_FLAGS "-march=armv8-a" CACHE STRING "Portable ARM64 architecture flags")

    message(STATUS "Detected ARM64 architecture (${CMAKE_SYSTEM_PROCESSOR})")
    message(STATUS "ARM NEON support: YES (guaranteed by ARMv8-A)")

elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64|x64)")
    set(ARCH_X86_64 ON CACHE BOOL "x86_64 architecture detected")
    set(ARCH_ARM64 OFF CACHE BOOL "ARM64 architecture not detected")

    # x86-64 guarantees SSE2 support
    add_compile_definitions(ARCH_X86_64=1)

    # Portable x86_64 flags
    set(ARCH_FLAGS "-msse2 -march=x86-64" CACHE STRING "Portable x86_64 architecture flags")

    message(STATUS "Detected x86_64 architecture (${CMAKE_SYSTEM_PROCESSOR})")
    message(STATUS "SSE2 support: YES (guaranteed by x86-64)")

else()
    set(ARCH_X86_64 OFF CACHE BOOL "x86_64 architecture not detected")
    set(ARCH_ARM64 OFF CACHE BOOL "ARM64 architecture not detected")

    message(WARNING "Unknown architecture: ${CMAKE_SYSTEM_PROCESSOR}")
    message(WARNING "Building without architecture-specific optimizations")

    set(ARCH_FLAGS "" CACHE STRING "No architecture-specific flags")
endif()

# Export variables for use in other CMake files
mark_as_advanced(ARCH_X86_64 ARCH_ARM64 ARCH_FLAGS)

