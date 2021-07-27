macro(get_sha1 args1 args2)
  execute_process(
    COMMAND git describe --always 
    WORKING_DIRECTORY ${args1}
    OUTPUT_VARIABLE ${args2}
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()