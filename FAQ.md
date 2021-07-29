## My machine with built-in gcc <6. Even if I set LD_LIBRARY_PATH and PATH to use new gcc, the source code still tried to use old gcc in /usr/bin and failed compilation.
The key is setting LD_LIBRARY_PATH and PATH before typing `cmake ./`. The easiest way maybe:  
1) Delete the current source code using commands `rm -rf anchorwave`.  
2) Set `LD_LIBRARY_PATH`, `PATH`, `CC`, and `GCC` to the new gcc.  
3) Re-clone the AnchorWave repository and compile it.  
## Should I perform genome repeat masking before feed into AnchorWave
Genome masking is not expected to improve the performance of AnchorWave.  
AnchorWave do not utilize any soft masking information. Hard masking would increase the
computational cost of AnchorWave.
