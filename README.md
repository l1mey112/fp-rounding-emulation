# https://stackoverflow.com/questions/78776730/double-floating-point-operations-using-four-ieee-rounding-modes-implemented-in-t/78779343#comment138931233_78779343

```
l-m@debdinger ~/g/fp-rounding-emulation (master)> gcc *.c -lm -Wshadow -O3 && ./a.out
fadd:
  fprc(0): 2clks
  fprc(1): 45clks [100000/100000] x17 overhead
  fprc(2): 21clks [100000/100000] x8 overhead
  fprc(3): 10clks [100000/100000] x4 overhead
fsub:
  fprc(0): 2clks
  fprc(1): 46clks [100000/100000] x17 overhead
  fprc(2): 22clks [100000/100000] x8 overhead
  fprc(3): 10clks [100000/100000] x4 overhead
fmul:
  fprc(0): 3clks
  fprc(1): 76clks [100000/100000] x23 overhead
  fprc(2): 38clks [100000/100000] x12 overhead
  fprc(3): 39clks [100000/100000] x12 overhead
fmul (fma):
  fprc(0): 3clks
  fprc(1): 32clks [100000/100000] x10 overhead
  fprc(2): 12clks [100000/100000] x4 overhead
  fprc(3): 14clks [100000/100000] x4 overhead
fdiv:
  fprc(0): 3clks
  fprc(1): 94clks [100000/100000] x24 overhead
  fprc(2): 64clks [100000/100000] x16 overhead
  fprc(3): 64clks [100000/100000] x16 overhead
fdiv (fma):
  fprc(0): 4clks
  fprc(1): 45clks [100000/100000] x10 overhead
  fprc(2): 18clks [100000/100000] x4 overhead
  fprc(3): 17clks [100000/100000] x4 overhead
fsqrt:
  fprc(0): 8clks
  fprc(1): 97clks [100000/100000] x12 overhead
  fprc(2): 67clks [100000/100000] x8 overhead
  fprc(3): 69clks [100000/100000] x8 overhead
fsqrt (fma):
  fprc(0): 8clks
  fprc(1): 52clks [100000/100000] x6 overhead
  fprc(2): 24clks [100000/100000] x3 overhead
  fprc(3): 23clks [100000/100000] x3 overhead
```