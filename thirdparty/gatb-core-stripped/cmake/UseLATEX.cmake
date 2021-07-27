# CMake utility to compile Latex documents with pdflatex 
# Version: 1.0.1
# Author: Baptiste Wicht <baptiste.wicht@gmail.com>

# Original statement
# Author: Kenneth Moreland <kmorel@sandia.gov>
#
# Copyright 2004 Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the
# U.S. Government. Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that this Notice and any
# statement of authorship are reproduced on all copies.

# The following function is defined:
# ADD_LATEX_DOCUMENT(<tex_file>
#                    [BIBFILES <bib_files>]
#                    [INPUTS <input_tex_files>]
#                    [IMAGE_DIRS] <image_directories>
#                    [IMAGES] <image_files>
#                    [CONFIGURE] <tex_files>
#                    [DEPENDS] <tex_files>
#                    [FILTER_OUTPUT]
#                    [USE_INDEX] 
#                    [USE_GLOSSARY])
#
# Adds targets that compile <tex_file>.The latex output is placed in LATEX_OUTPUT_PATH 
# or CMAKE_CURRENT_BINARY_DIR if the former is not set.

# Change log
#
# Version 1.0.1
#
# Add filter feature
#
# Version 1.0.0
#
# Clean up version of Kenneth Moreland

#############################################################################
# Find the location of myself while originally executing.  If you do this
# inside of a macro, it will recode where the macro was invoked.
#############################################################################
SET(LATEX_USE_LATEX_LOCATION ${CMAKE_CURRENT_LIST_FILE}
    CACHE INTERNAL "Location of UseLATEX.cmake file." FORCE
    )

#############################################################################
# Generic helper functions
#############################################################################

FUNCTION(LATEX_LIST_CONTAINS var value)
    SET(input_list ${ARGN})
    LIST(FIND input_list "${value}" index)
    IF (index GREATER -1)
        SET(${var} TRUE PARENT_SCOPE)
    ELSE (index GREATER -1)
        SET(${var} PARENT_SCOPE)
    ENDIF (index GREATER -1)
ENDFUNCTION(LATEX_LIST_CONTAINS)

# Parse function arguments.  Variables containing the results are placed
# in the global scope for historical reasons.
FUNCTION(LATEX_PARSE_ARGUMENTS prefix arg_names option_names)
    SET(DEFAULT_ARGS)
    FOREACH(arg_name ${arg_names})
        SET(${prefix}_${arg_name} CACHE INTERNAL "${prefix} argument" FORCE)
    ENDFOREACH(arg_name)
    FOREACH(option ${option_names})
        SET(${prefix}_${option} CACHE INTERNAL "${prefix} option" FORCE)
    ENDFOREACH(option)

    SET(current_arg_name DEFAULT_ARGS)
    SET(current_arg_list)
    FOREACH(arg ${ARGN})
        LATEX_LIST_CONTAINS(is_arg_name ${arg} ${arg_names})
        LATEX_LIST_CONTAINS(is_option ${arg} ${option_names})
        IF (is_arg_name)
            SET(${prefix}_${current_arg_name} ${current_arg_list}
                CACHE INTERNAL "${prefix} argument" FORCE)
            SET(current_arg_name ${arg})
            SET(current_arg_list)
        ELSEIF (is_option)
            SET(${prefix}_${arg} TRUE CACHE INTERNAL "${prefix} option" FORCE)
        ELSE (is_arg_name)
            SET(current_arg_list ${current_arg_list} ${arg})
        ENDIF (is_arg_name)
    ENDFOREACH(arg)
    SET(${prefix}_${current_arg_name} ${current_arg_list}
        CACHE INTERNAL "${prefix} argument" FORCE)
ENDFUNCTION(LATEX_PARSE_ARGUMENTS)

# Match the contents of a file to a regular expression.
FUNCTION(LATEX_FILE_MATCH variable filename regexp default)
    # The FILE STRINGS command would be a bit better, but I'm not totally sure
    # the match will always be to a whole line, and I don't want to break things.
    FILE(READ ${filename} file_contents)
    STRING(REGEX MATCHALL "${regexp}"
        match_result ${file_contents}
        )
    IF (match_result)
        SET(${variable} "${match_result}" PARENT_SCOPE)
    ELSE (match_result)
        SET(${variable} "${default}" PARENT_SCOPE)
    ENDIF (match_result)
ENDFUNCTION(LATEX_FILE_MATCH)

# A version of GET_FILENAME_COMPONENT that treats extensions after the last
# period rather than the first.  To the best of my knowledge, all filenames
# typically used by LaTeX, including image files, have small extensions
# after the last dot.
FUNCTION(LATEX_GET_FILENAME_COMPONENT varname filename type)
    SET(result)
    IF ("${type}" STREQUAL "NAME_WE")
        GET_FILENAME_COMPONENT(name ${filename} NAME)
        STRING(REGEX REPLACE "\\.[^.]*\$" "" result "${name}")
    ELSEIF ("${type}" STREQUAL "EXT")
        GET_FILENAME_COMPONENT(name ${filename} NAME)
        STRING(REGEX MATCH "\\.[^.]*\$" result "${name}")
    ELSE ("${type}" STREQUAL "NAME_WE")
        GET_FILENAME_COMPONENT(result ${filename} ${type})
    ENDIF ("${type}" STREQUAL "NAME_WE")
    SET(${varname} "${result}" PARENT_SCOPE)
ENDFUNCTION(LATEX_GET_FILENAME_COMPONENT)

#############################################################################
# Functions that perform processing during a LaTeX build.
#############################################################################
FUNCTION(LATEX_MAKEGLOSSARIES)
    # This is really a bare bones port of the makeglossaries perl script into
    # CMake scripting.
    IF (NOT LATEX_TARGET)
        MESSAGE(SEND_ERROR "Need to define LATEX_TARGET")
    ENDIF (NOT LATEX_TARGET)

    SET(aux_file ${LATEX_TARGET}.aux)

    IF (NOT EXISTS ${aux_file})
        MESSAGE(SEND_ERROR "${aux_file} does not exist.  Run latex on your target file.")
    ENDIF (NOT EXISTS ${aux_file})

    LATEX_FILE_MATCH(newglossary_lines ${aux_file}
        "@newglossary[ \t]*{([^}]*)}{([^}]*)}{([^}]*)}{([^}]*)}"
        "@newglossary{main}{glg}{gls}{glo}"
        )

    LATEX_FILE_MATCH(istfile_line ${aux_file}
        "@istfilename[ \t]*{([^}]*)}"
        "@istfilename{${LATEX_TARGET}.ist}"
        )
    STRING(REGEX REPLACE "@istfilename[ \t]*{([^}]*)}" "\\1"
        istfile ${istfile_line}
        )

    IF (NOT MAKEINDEX_COMPILER)
        MESSAGE(SEND_ERROR "Need to define MAKEINDEX_COMPILER")
    ENDIF (NOT MAKEINDEX_COMPILER)

    FOREACH(newglossary ${newglossary_lines})
        STRING(REGEX REPLACE
            "@newglossary[ \t]*{([^}]*)}{([^}]*)}{([^}]*)}{([^}]*)}"
            "\\1" glossary_name ${newglossary}
            )
        STRING(REGEX REPLACE
            "@newglossary[ \t]*{([^}]*)}{([^}]*)}{([^}]*)}{([^}]*)}"
            "${LATEX_TARGET}.\\2" glossary_log ${newglossary}
            )
        STRING(REGEX REPLACE
            "@newglossary[ \t]*{([^}]*)}{([^}]*)}{([^}]*)}{([^}]*)}"
            "${LATEX_TARGET}.\\3" glossary_out ${newglossary}
            )
        STRING(REGEX REPLACE
            "@newglossary[ \t]*{([^}]*)}{([^}]*)}{([^}]*)}{([^}]*)}"
            "${LATEX_TARGET}.\\4" glossary_in ${newglossary}
            )

        IF (LATEX_FILTER_OUTPUT)
            EXECUTE_PROCESS(
                COMMAND ${MAKEINDEX_COMPILER} ${MAKEGLOSSARIES_COMPILER_FLAGS} -s ${istfile} -t ${glossary_log} -o ${glossary_out} ${glossary_in}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/.tmp.log
                ERROR_FILE ${CMAKE_CURRENT_SOURCE_DIR}/.tmp.log
                )
            EXECUTE_PROCESS(
                COMMAND awk -f ../index_filter.awk ${CMAKE_CURRENT_SOURCE_DIR}/.tmp.log
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                )
        ELSE ()
            EXECUTE_PROCESS(
                COMMAND ${MAKEINDEX_COMPILER} ${MAKEGLOSSARIES_COMPILER_FLAGS} -s ${istfile} -t ${glossary_log} -o ${glossary_out} ${glossary_in}
                WORKING_DIRECTORY ${LATEX_OUTPUT}
                )
        ENDIF (LATEX_FILTER_OUTPUT)

    ENDFOREACH(newglossary)
ENDFUNCTION(LATEX_MAKEGLOSSARIES)

#############################################################################
# Helper functions for establishing LaTeX build.
#############################################################################

FUNCTION(LATEX_NEEDIT VAR NAME)
    IF (NOT ${VAR})
        MESSAGE(SEND_ERROR "I need the ${NAME} command.")
    ENDIF(NOT ${VAR})
ENDFUNCTION(LATEX_NEEDIT)

FUNCTION(LATEX_SETUP_VARIABLES)
    SET(LATEX_OUTPUT_PATH "${LATEX_OUTPUT_PATH}"
        CACHE PATH "If non empty, specifies the location to place LaTeX output."
        )

    FIND_PACKAGE(LATEX)
    FIND_PACKAGE(UnixCommands)

    MARK_AS_ADVANCED(CLEAR
        LATEX_COMPILER
        PDFLATEX_COMPILER
        BIBTEX_COMPILER
        MAKEINDEX_COMPILER
        )

    LATEX_NEEDIT(PDFLATEX_COMPILER pdflatex)
    LATEX_NEEDIT(BIBTEX_COMPILER bibtex)
    LATEX_NEEDIT(MAKEINDEX_COMPILER makeindex)

    SET(LATEX_COMPILER_FLAGS "-interaction=errorstopmode"
        CACHE STRING "Flags passed to latex.")
    SET(PDFLATEX_COMPILER_FLAGS ${LATEX_COMPILER_FLAGS}
        CACHE STRING "Flags passed to pdflatex.")
    SET(BIBTEX_COMPILER_FLAGS ""
        CACHE STRING "Flags passed to bibtex.")
    SET(MAKEINDEX_COMPILER_FLAGS ""
        CACHE STRING "Flags passed to makeindex.")
    SET(MAKEGLOSSARIES_COMPILER_FLAGS ""
        CACHE STRING "Flags passed to makeglossaries.")
    MARK_AS_ADVANCED(
        LATEX_COMPILER_FLAGS
        PDFLATEX_COMPILER_FLAGS
        BIBTEX_COMPILER_FLAGS
        MAKEINDEX_COMPILER_FLAGS
        MAKEGLOSSARIES_COMPILER_FLAGS
        )
    SEPARATE_ARGUMENTS(LATEX_COMPILER_FLAGS)
    SEPARATE_ARGUMENTS(PDFLATEX_COMPILER_FLAGS)
    SEPARATE_ARGUMENTS(BIBTEX_COMPILER_FLAGS)
    SEPARATE_ARGUMENTS(MAKEINDEX_COMPILER_FLAGS)
    SEPARATE_ARGUMENTS(MAKEGLOSSARIES_COMPILER_FLAGS)

    FIND_PROGRAM(IMAGEMAGICK_CONVERT convert
        DOC "The convert program that comes with ImageMagick (available at http://www.imagemagick.org)."
        )

    IF (NOT IMAGEMAGICK_CONVERT)
        MESSAGE(SEND_ERROR "Could not find convert program.  Please download ImageMagick from http://www.imagemagick.org and install.")
    ENDIF (NOT IMAGEMAGICK_CONVERT)
    
    FIND_PROGRAM(CAIRO_CONVERT cairosvg
        DOC "Cairo SVG Converter"
        )

    SET(LATEX_RASTER_SCALE 100)
    SET(LATEX_OPPOSITE_RASTER_SCALE 16)

    SET(LATEX_SVG_IMAGE_EXTENSIONS .svg CACHE INTERNAL "")
    SET(LATEX_EPS_IMAGE_EXTENSIONS .eps CACHE INTERNAL "")
    SET(LATEX_PDF_VECTOR_IMAGE_EXTENSIONS .pdf CACHE INTERNAL "")
    SET(LATEX_PDF_RASTER_IMAGE_EXTENSIONS .png .jpeg .jpg CACHE INTERNAL "")
    SET(LATEX_PDF_IMAGE_EXTENSIONS ${LATEX_PDF_VECTOR_IMAGE_EXTENSIONS} ${LATEX_PDF_RASTER_IMAGE_EXTENSIONS} CACHE INTERNAL "")
    SET(LATEX_IMAGE_EXTENSIONS ${LATEX_PDF_IMAGE_EXTENSIONS} ${LATEX_EPS_IMAGE_EXTENSIONS} ${LATEX_SVG_IMAGE_EXTENSIONS} CACHE INTERNAL "")
ENDFUNCTION(LATEX_SETUP_VARIABLES)

FUNCTION(LATEX_GET_OUTPUT_PATH var)
    SET(latex_output_path)
    IF (LATEX_OUTPUT_PATH)
        IF ("${LATEX_OUTPUT_PATH}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
            MESSAGE(SEND_ERROR "You cannot set LATEX_OUTPUT_PATH to the same directory that contains LaTeX input files.")
        ELSE ("${LATEX_OUTPUT_PATH}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
            SET(latex_output_path "${LATEX_OUTPUT_PATH}")
        ENDIF ("${LATEX_OUTPUT_PATH}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    ELSE (LATEX_OUTPUT_PATH)
        IF ("${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
            MESSAGE(SEND_ERROR "LaTeX files must be built out of source or you must set LATEX_OUTPUT_PATH.")
        ELSE ("${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
            SET(latex_output_path "${CMAKE_CURRENT_BINARY_DIR}")
        ENDIF ("${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    ENDIF (LATEX_OUTPUT_PATH)
    SET(${var} ${latex_output_path} PARENT_SCOPE)
ENDFUNCTION(LATEX_GET_OUTPUT_PATH)

FUNCTION(LATEX_ADD_CONVERT_COMMAND
        output_path
        input_path
        output_extension
        input_extension
        flags
        )

    IF (${input_extension} STREQUAL ".svg" AND ${output_extension} STREQUAL ".pdf")
        ADD_CUSTOM_COMMAND(OUTPUT ${output_path}
            COMMAND ${CAIRO_CONVERT}
            ARGS ${input_path} "-o" ${output_path}
            DEPENDS ${input_path}
            )
    ELSE ()
        SET (converter ${IMAGEMAGICK_CONVERT})
        SET (convert_flags "")

        IF (${input_extension} STREQUAL ".eps" AND ${output_extension} STREQUAL ".pdf")
            # ImageMagick has broken eps to pdf conversion, use ps2pdf instead
            IF (PS2PDF_CONVERTER)
                SET (converter ${PS2PDF_CONVERTER})
                SET (convert_flags -dEPSCrop ${PS2PDF_CONVERTER_FLAGS})
            ELSE (PS2PDF_CONVERTER)
                MESSAGE(SEND_ERROR "Using postscript files with pdflatex requires ps2pdf for conversion.")
            ENDIF (PS2PDF_CONVERTER)
        ELSE (${input_extension} STREQUAL ".eps" AND ${output_extension} STREQUAL ".pdf")
            SET (convert_flags ${flags})
        ENDIF (${input_extension} STREQUAL ".eps" AND ${output_extension} STREQUAL ".pdf")

        ADD_CUSTOM_COMMAND(OUTPUT ${output_path}
            COMMAND ${converter}
            ARGS ${convert_flags} ${input_path} ${output_path}
            DEPENDS ${input_path}
            )
    ENDIF()
ENDFUNCTION(LATEX_ADD_CONVERT_COMMAND)

# Makes custom commands to convert a file to a particular type.
FUNCTION(LATEX_CONVERT_IMAGE
        output_files_var
        input_file
        output_extension
        convert_flags
        output_extensions
        other_files
        )
    
    SET(output_file_list)
    SET(input_dir ${CMAKE_CURRENT_SOURCE_DIR})
    
    LATEX_GET_OUTPUT_PATH(output_dir)
    LATEX_GET_FILENAME_COMPONENT(extension "${input_file}" EXT)

    # Check input filename for potential problems with LaTeX.
    LATEX_GET_FILENAME_COMPONENT(name "${input_file}" NAME_WE)
    IF (name MATCHES ".*\\..*")
        STRING(REPLACE "." "-" suggested_name "${name}")
        SET(suggested_name "${suggested_name}${extension}")
        MESSAGE(WARNING "Some LaTeX distributions have problems with image file names with multiple extensions.  Consider changing ${name}${extension} to something like ${suggested_name}.")
    ENDIF (name MATCHES ".*\\..*")

    STRING(REGEX REPLACE "\\.[^.]*\$" ${output_extension} output_file "${input_file}")

    LATEX_LIST_CONTAINS(is_type ${extension} ${output_extensions})
    
    IF (is_type)
        ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/${input_file}
            COMMAND ${CMAKE_COMMAND}
            ARGS -E copy ${input_dir}/${input_file} ${output_dir}/${input_file}
            DEPENDS ${input_dir}/${input_file}
            )
        SET(output_file_list ${output_file_list} ${output_dir}/${input_file})
    ELSE (is_type)
        LATEX_ADD_CONVERT_COMMAND(${output_dir}/${output_file}
            ${input_dir}/${input_file} ${output_extension} ${extension}
            "${convert_flags}")
        SET(output_file_list ${output_file_list} ${output_dir}/${output_file})
    ENDIF (is_type)

    SET(${output_files_var} ${output_file_list} PARENT_SCOPE)
ENDFUNCTION(LATEX_CONVERT_IMAGE)

# Adds custom commands to process the given files for pdf builds.
# Adds the output files to the given variables (does not replace).
FUNCTION(LATEX_PROCESS_IMAGES pdf_outputs_var)
    LATEX_GET_OUTPUT_PATH(output_dir)
    SET(pdf_outputs)
    FOREACH(file ${ARGN})
        IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
            LATEX_GET_FILENAME_COMPONENT(extension "${file}" EXT)
            SET(convert_flags)

            LATEX_LIST_CONTAINS(is_raster "${extension}" ${LATEX_PDF_RASTER_IMAGE_EXTENSIONS})
            LATEX_LIST_CONTAINS(is_svg "${extension}" ${LATEX_SVG_IMAGE_EXTENSIONS})

            # Make sure the output directory exists.
            LATEX_GET_FILENAME_COMPONENT(path "${output_dir}/${file}" PATH)
            MAKE_DIRECTORY("${path}")

            # Do conversions for PDF and SVG
            IF (is_raster)
                LATEX_CONVERT_IMAGE(output_files "${file}" .png "${convert_flags}"
                    "${LATEX_PDF_IMAGE_EXTENSIONS}" "${ARGN}")
                SET(pdf_outputs ${pdf_outputs} ${output_files})
            ELSEIF (is_svg)
                LATEX_CONVERT_IMAGE(output_files "${file}" .pdf "${convert_flags}"
                    "${LATEX_PDF_IMAGE_EXTENSIONS}" "${ARGN}")
                SET(pdf_outputs ${pdf_outputs} ${output_files})
            ELSE ()
                LATEX_CONVERT_IMAGE(output_files "${file}" .pdf "${convert_flags}"
                    "${LATEX_PDF_IMAGE_EXTENSIONS}" "${ARGN}")
                SET(pdf_outputs ${pdf_outputs} ${output_files})
            ENDIF ()
        ELSE (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
            MESSAGE(WARNING "Could not find file ${CMAKE_CURRENT_SOURCE_DIR}/${file}.  Are you sure you gave relative paths to IMAGES?")
        ENDIF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
    ENDFOREACH(file)

    SET(${pdf_outputs_var} ${pdf_outputs} PARENT_SCOPE)
ENDFUNCTION(LATEX_PROCESS_IMAGES)

FUNCTION(LATEX_COPY_GLOBBED_FILES pattern dest)
    FILE(GLOB file_list ${pattern})
    FOREACH(in_file ${file_list})
        LATEX_GET_FILENAME_COMPONENT(out_file ${in_file} NAME)
        CONFIGURE_FILE(${in_file} ${dest}/${out_file} COPYONLY)
    ENDFOREACH(in_file)
ENDFUNCTION(LATEX_COPY_GLOBBED_FILES)

FUNCTION(LATEX_COPY_INPUT_FILE file)
    LATEX_GET_OUTPUT_PATH(output_dir)

    IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${file})
        LATEX_GET_FILENAME_COMPONENT(path ${file} PATH)
        FILE(MAKE_DIRECTORY ${output_dir}/${path})

        LATEX_LIST_CONTAINS(use_config ${file} ${LATEX_CONFIGURE})
        IF (use_config)
            CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${file}
                ${output_dir}/${file}
                @ONLY
                )
            ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/${file}
                COMMAND ${CMAKE_COMMAND}
                ARGS ${CMAKE_BINARY_DIR}
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${file}
                )
        ELSE (use_config)
            ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/${file}
                COMMAND ${CMAKE_COMMAND}
                ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${output_dir}/${file}
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${file}
                )
        ENDIF (use_config)
    ELSE (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${file})
        IF (EXISTS ${output_dir}/${file})
            # Special case: output exists but input does not.  Assume that it was
            # created elsewhere and skip the input file copy.
        ELSE (EXISTS ${output_dir}/${file})
            MESSAGE("Could not find input file ${CMAKE_CURRENT_SOURCE_DIR}/${file}")
        ENDIF (EXISTS ${output_dir}/${file})
    ENDIF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${file})
ENDFUNCTION(LATEX_COPY_INPUT_FILE)

#############################################################################
# Commands provided by the UseLATEX.cmake "package"
#############################################################################

FUNCTION(LATEX_USAGE command message)
    MESSAGE(SEND_ERROR
        "${message}\nUsage: ${command}(<tex_file>\n           [BIBFILES <bib_file> <bib_file> ...]\n           [INPUTS <tex_file> <tex_file> ...]\n           [IMAGE_DIRS <directory1> <directory2> ...]\n           [IMAGES <image_file1> <image_file2>\n           [CONFIGURE <tex_file> <tex_file> ...]\n           [DEPENDS <tex_file> <tex_file> ...]\n           [USE_INDEX] [FILTER_OUTPUT] [USE_GLOSSARY])"
        )
ENDFUNCTION(LATEX_USAGE command message)

# Parses arguments to ADD_LATEX_DOCUMENT and ADD_LATEX_TARGETS and sets the
# variables LATEX_TARGET, LATEX_IMAGE_DIR, LATEX_BIBFILES, LATEX_DEPENDS, and
# LATEX_INPUTS.
FUNCTION(PARSE_ADD_LATEX_ARGUMENTS command)
    LATEX_PARSE_ARGUMENTS(
        LATEX
        "BIBFILES;INPUTS;IMAGE_DIRS;IMAGES;CONFIGURE;DEPENDS"
        "USE_INDEX;FILTER_OUTPUT;USE_GLOSSARY;USE_GLOSSARIES"
        ${ARGN}
        )

    # The first argument is the target latex file.
    IF (LATEX_DEFAULT_ARGS)
        LIST(GET LATEX_DEFAULT_ARGS 0 latex_main_input)
        LIST(REMOVE_AT LATEX_DEFAULT_ARGS 0)
        LATEX_GET_FILENAME_COMPONENT(latex_target ${latex_main_input} NAME_WE)
        SET(LATEX_MAIN_INPUT ${latex_main_input} CACHE INTERNAL "" FORCE)
        SET(LATEX_TARGET ${latex_target} CACHE INTERNAL "" FORCE)
    ELSE (LATEX_DEFAULT_ARGS)
        LATEX_USAGE(${command} "No tex file target given to ${command}.")
    ENDIF (LATEX_DEFAULT_ARGS)

    IF (LATEX_DEFAULT_ARGS)
        LATEX_USAGE(${command} "Invalid or depricated arguments: ${LATEX_DEFAULT_ARGS}")
    ENDIF (LATEX_DEFAULT_ARGS)

    # Backward compatibility between 1.6.0 and 1.6.1.
    IF (LATEX_USE_GLOSSARIES)
        SET(LATEX_USE_GLOSSARY TRUE CACHE INTERNAL "" FORCE)
    ENDIF (LATEX_USE_GLOSSARIES)
ENDFUNCTION(PARSE_ADD_LATEX_ARGUMENTS)

FUNCTION(ADD_LATEX_TARGETS_INTERNAL)
    # The commands to run LaTeX.  They are repeated multiple times.
    SET(pdflatex_draft_command ${PDFLATEX_COMPILER} -draftmode -shell-escape ${PDFLATEX_COMPILER_FLAGS} ${LATEX_MAIN_INPUT})
    SET(pdflatex_build_command ${PDFLATEX_COMPILER} -shell-escape ${PDFLATEX_COMPILER_FLAGS} ${LATEX_MAIN_INPUT})

    # The command to create the index
    SET(makeindex_command ${MAKEINDEX_COMPILER} ${MAKEINDEX_COMPILER_FLAGS} ${LATEX_TARGET}.idx)

    # The command to create the bibliography
    SET(bibtex_command ${BIBTEX_COMPILER} ${BIBTEX_COMPILER_FLAGS} ${LATEX_TARGET})

    IF (LATEX_FILTER_OUTPUT)                                       
        SET(pdflatex_draft_command ${pdflatex_draft_command} | awk -f reverse.awk | awk -f compose.awk | awk -f reverse.awk | sed \"s/\\\\[.\\+\\\\]//\" | awk -f filter.awk)
        SET(pdflatex_build_command ${pdflatex_build_command} | awk -f reverse.awk | awk -f compose.awk | awk -f reverse.awk | sed \"s/\\\\[.\\+\\\\]//\" | awk -f filter.awk)
        
        SET(makeindex_command ${makeindex_command} | awk -f index_filter.awk)
        SET(bibtex_command ${bibtex_command} | awk -f bibtex_filter.awk)
    ENDIF (LATEX_FILTER_OUTPUT)

    # Set up target names.
    SET(pdf_target      pdf)
    SET(auxclean_target auxclean)

    # Probably not all of these will be generated, but they could be.
    # Note that the aux file is added later.
    SET(auxiliary_clean_files
        ${output_dir}/${LATEX_TARGET}.bbl
        ${output_dir}/${LATEX_TARGET}.blg
        ${output_dir}/${LATEX_TARGET}-blx.bib
        ${output_dir}/${LATEX_TARGET}.glg
        ${output_dir}/${LATEX_TARGET}.glo
        ${output_dir}/${LATEX_TARGET}.gls
        ${output_dir}/${LATEX_TARGET}.idx
        ${output_dir}/${LATEX_TARGET}.ilg
        ${output_dir}/${LATEX_TARGET}.ind
        ${output_dir}/${LATEX_TARGET}.ist
        ${output_dir}/${LATEX_TARGET}.log
        ${output_dir}/${LATEX_TARGET}.lol
        ${output_dir}/${LATEX_TARGET}.tdo
        ${output_dir}/${LATEX_TARGET}.out
        ${output_dir}/${LATEX_TARGET}.toc
        ${output_dir}/${LATEX_TARGET}.lof
        ${output_dir}/${LATEX_TARGET}.xdy
        ${output_dir}/${LATEX_TARGET}.dvi
        ${output_dir}/${LATEX_TARGET}.ps
        ${output_dir}/${LATEX_TARGET}.pdf
        )

    SET(image_list ${LATEX_IMAGES})

    # For each directory in LATEX_IMAGE_DIRS, glob all the image files and
    # place them in LATEX_IMAGES.
    FOREACH(dir ${LATEX_IMAGE_DIRS})
        IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
            MESSAGE(WARNING "Image directory ${CMAKE_CURRENT_SOURCE_DIR}/${dir} does not exist.  Are you sure you gave relative directories to IMAGE_DIRS?")
        ENDIF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
        
        FOREACH(extension ${LATEX_IMAGE_EXTENSIONS})
            FILE(GLOB files ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*${extension})
            FOREACH(file ${files})
                LATEX_GET_FILENAME_COMPONENT(filename ${file} NAME)
                SET(image_list ${image_list} ${dir}/${filename})
            ENDFOREACH(file)
        ENDFOREACH(extension)
    ENDFOREACH(dir)

    LATEX_PROCESS_IMAGES(pdf_images ${image_list})

    SET(make_pdf_command 
        ${CMAKE_COMMAND} -E chdir ${output_dir}
        ${pdflatex_draft_command}
        )

    SET(make_pdf_depends ${LATEX_DEPENDS} ${pdf_images})
    FOREACH(input ${LATEX_MAIN_INPUT} ${LATEX_INPUTS})
        SET(make_pdf_depends ${make_pdf_depends} ${output_dir}/${input})
        IF (${input} MATCHES "\\.tex$")
            STRING(REGEX REPLACE "\\.tex$" "" input_we ${input})
            SET(auxiliary_clean_files ${auxiliary_clean_files}
                ${output_dir}/${input_we}.aux
                ${output_dir}/${input}.aux
                )
        ENDIF (${input} MATCHES "\\.tex$")
    ENDFOREACH(input)

    IF (LATEX_USE_GLOSSARY)
        FOREACH(dummy 0 1)   # Repeat these commands twice.
            SET(make_pdf_command ${make_pdf_command}
                COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
                ${CMAKE_COMMAND}
                -D LATEX_BUILD_COMMAND=makeglossaries
                -D LATEX_TARGET=${LATEX_TARGET}
                -D MAKEINDEX_COMPILER=${MAKEINDEX_COMPILER}
                -D MAKEGLOSSARIES_COMPILER_FLAGS=${MAKEGLOSSARIES_COMPILER_FLAGS}
                -D LATEX_OUTPUT=${output_dir}
                -D LATEX_FILTER_OUTPUT=${LATEX_FILTER_OUTPUT}
                -P ${LATEX_USE_LATEX_LOCATION}
                COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
                ${pdflatex_draft_command}
                )
        ENDFOREACH(dummy)
    ENDIF (LATEX_USE_GLOSSARY)

    IF (LATEX_BIBFILES)
        SET(make_pdf_command ${make_pdf_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${bibtex_command}
            )

        FOREACH (bibfile ${LATEX_BIBFILES})
            SET(make_pdf_depends ${make_pdf_depends} ${output_dir}/${bibfile})
        ENDFOREACH (bibfile ${LATEX_BIBFILES})
    ENDIF (LATEX_BIBFILES)

    IF (LATEX_USE_INDEX)
        SET(make_pdf_command ${make_pdf_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir} ${pdflatex_draft_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir} ${makeindex_command} 
            )
    ENDIF (LATEX_USE_INDEX)

    # In fast mode, only do one pass
    SET(make_pdf_fast_command ${make_pdf_fast_command}
        COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
        ${pdflatex_build_command})
    
    # Do two pass at the end
    SET(make_pdf_command ${make_pdf_command}
        COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
        ${pdflatex_draft_command}
        COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
        ${pdflatex_build_command})

    # If all features are enabled, create faster version
    IF (LATEX_USE_INDEX AND LATEX_BIBFILES AND LATEX_USE_GLOSSARY)
        SET(make_pdf_command
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${pdflatex_draft_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${bibtex_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${pdflatex_draft_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${makeindex_command} 
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${CMAKE_COMMAND}
            -D LATEX_BUILD_COMMAND=makeglossaries
            -D LATEX_TARGET=${LATEX_TARGET}
            -D MAKEINDEX_COMPILER=${MAKEINDEX_COMPILER}
            -D MAKEGLOSSARIES_COMPILER_FLAGS=${MAKEGLOSSARIES_COMPILER_FLAGS}
            -D LATEX_OUTPUT=${output_dir}
            -D LATEX_FILTER_OUTPUT=${LATEX_FILTER_OUTPUT}
            -P ${LATEX_USE_LATEX_LOCATION}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${pdflatex_draft_command}
            COMMAND ${CMAKE_COMMAND} -E chdir ${output_dir}
            ${pdflatex_build_command})
    ENDIF (LATEX_USE_INDEX AND LATEX_BIBFILES AND LATEX_USE_GLOSSARY)

    # Finally add the target to the makefile
    ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/${LATEX_TARGET}.pdf
        COMMAND ${make_pdf_command}
        DEPENDS ${make_pdf_depends}
        )
    
    ADD_CUSTOM_TARGET(${pdf_target} ALL
        DEPENDS ${output_dir}/${LATEX_TARGET}.pdf)

    ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/fast_${LATEX_TARGET}.pdf
        COMMAND ${make_pdf_fast_command}
        DEPENDS ${make_pdf_depends}
        )
    ADD_CUSTOM_TARGET(fast
        DEPENDS ${output_dir}/fast_${LATEX_TARGET}.pdf)

    SET_DIRECTORY_PROPERTIES(.
        ADDITIONAL_MAKE_CLEAN_FILES "${auxiliary_clean_files}"
        )

    ADD_CUSTOM_TARGET(${auxclean_target}
        COMMENT "Cleaning auxiliary LaTeX files."
        COMMAND ${CMAKE_COMMAND} -E remove ${auxiliary_clean_files}
        )
ENDFUNCTION(ADD_LATEX_TARGETS_INTERNAL)

FUNCTION(ADD_LATEX_TARGETS)
    LATEX_GET_OUTPUT_PATH(output_dir)
    PARSE_ADD_LATEX_ARGUMENTS(ADD_LATEX_TARGETS ${ARGV})

    ADD_LATEX_TARGETS_INTERNAL()
ENDFUNCTION(ADD_LATEX_TARGETS)

FUNCTION(ADD_LATEX_DOCUMENT)
    LATEX_GET_OUTPUT_PATH(output_dir)
    IF (output_dir)
        PARSE_ADD_LATEX_ARGUMENTS(ADD_LATEX_DOCUMENT ${ARGV})

        LATEX_COPY_INPUT_FILE(${LATEX_MAIN_INPUT})

        FOREACH (bib_file ${LATEX_BIBFILES})
            LATEX_COPY_INPUT_FILE(${bib_file})
        ENDFOREACH (bib_file)

        FOREACH (input ${LATEX_INPUTS})
            LATEX_COPY_INPUT_FILE(${input})
        ENDFOREACH(input)

        LATEX_COPY_GLOBBED_FILES(${CMAKE_CURRENT_SOURCE_DIR}/*.cls ${output_dir})
        LATEX_COPY_GLOBBED_FILES(${CMAKE_CURRENT_SOURCE_DIR}/*.bst ${output_dir})
        LATEX_COPY_GLOBBED_FILES(${CMAKE_CURRENT_SOURCE_DIR}/*.clo ${output_dir})
        LATEX_COPY_GLOBBED_FILES(${CMAKE_CURRENT_SOURCE_DIR}/*.sty ${output_dir})
        LATEX_COPY_GLOBBED_FILES(${CMAKE_CURRENT_SOURCE_DIR}/*.ist ${output_dir})

        ADD_LATEX_TARGETS_INTERNAL()
    ENDIF (output_dir)
ENDFUNCTION(ADD_LATEX_DOCUMENT)

#############################################################################
# Actually do stuff
#############################################################################

IF (LATEX_BUILD_COMMAND)
    SET(command_handled)

    IF ("${LATEX_BUILD_COMMAND}" STREQUAL makeglossaries)
        LATEX_MAKEGLOSSARIES()
        SET(command_handled TRUE)
    ENDIF ("${LATEX_BUILD_COMMAND}" STREQUAL makeglossaries)

    IF (NOT command_handled)
        MESSAGE(SEND_ERROR "Unknown command: ${LATEX_BUILD_COMMAND}")
    ENDIF (NOT command_handled)

ELSE (LATEX_BUILD_COMMAND)
    # Must be part of the actual configure (included from CMakeLists.txt).
    LATEX_SETUP_VARIABLES()
ENDIF (LATEX_BUILD_COMMAND)
