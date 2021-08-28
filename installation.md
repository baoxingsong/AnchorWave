AnchorWave uses the SIMD instructions to seep up the dynamic programming algorithm. Specific functions have been implemented for SSE2, SSE4.1, AVX2 and AVX512 instruction sets.  
Generally speaking, the time cost: SSE2 >= SSE4.1 > AVX2 > AVX512.  
To check what CPU instructions are supported by your machine, you could run this command:
```cat /proc/cpuinfo | grep "flags" | uniq```  
By default, we assume the machine supports SSE4.1.

## If you are using machine supports AVX512 and would like to take the advantage of that
Clone the repository, and replace the default CMakeLists.txt with the avx512 one.
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
mv CMakeLists_avx512.txt CMakeLists.txt
cmake ./
make
```


## If you are using machine supports AVX2 and would like to take the advantage of AVX2
Clone the repository, and replace the default CMakeLists.txt with the avx2 one.
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
mv CMakeLists_avx2.txt CMakeLists.txt
cmake ./
make
```


## If you are using old machine supports SSE2 but not SSE4.1

Clone the repository, and replace the default CMakeLists.txt with the sse2 one.
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
mv CMakeLists_sse2.txt CMakeLists.txt
cmake ./
make
```

## If you are using very ole CPU or CPU with other architecture (i.e. ARM)
AnchorWave was not tested on that kind of platform. Highly would not work.

