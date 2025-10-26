Take a look at `benches.txt` to check the benches. To compile and check:

```sh
clang -ggdb -lm -march=native -O3 bench.c soft.c semi.c -o bench && ./bench
```

(you must have fma, check with `lscpu | grep fma`)