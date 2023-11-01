# xorshift

This little project of mine is a sort of total suit for xorshifts. It will provide tools capable of calculating coefficients for xorshifts of any size and period

## Quick start
compile and run the project with `make`:
```
> make
> mat.exe
```

Other compile targets for make are:
* `fast` : default, compiled with `-O2`
* `plain`
* `db` : debug info
* `pre` : only transform with C preprocessor

Which creates a program named `mat.exe`
