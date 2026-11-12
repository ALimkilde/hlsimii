# Writes the web page's webbing catalogue: the full catalogue up to its "# cli only" line (all of it if there is none).
#   cmake -DIN=<full webbings.yaml> -DOUT=<page's webbings.yaml> -P webbings_subset.cmake
file(READ "${IN}" text)
string(FIND "${text}" "\n# cli only" cut)
if(cut GREATER -1)
  string(SUBSTRING "${text}" 0 ${cut} text)
  string(APPEND text "\n")
endif()
file(WRITE "${OUT}" "${text}")
