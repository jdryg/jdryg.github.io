---
title: "jitcc microblog #2: Csmith"
date:  2026-04-22 12:00:00 +0200
categories: compilers
tags: debugging bugs
---

## Csmith

[Csmith](https://github.com/csmith-project/csmith) is a random generator of C programs. You execute it and it gives you back a C file. You build it with 2 different compilers, run both executables and check their outputs. If they match, everything is (probably) OK. If they don't, there must be a bug in one of the compilers.

In my case, the one of the 2 compilers is my own. The other is Clang 22.1.0 using the following options:

```c
const char* clangCmdLine = "clang.exe "
    "-O1 "
    "-fno-ms-compatibility "
    "-o %s %s";
```

where `%s` is the Csmith generated file.

For the purpose of this experiment, I assume Clang to be bug-free. So whatever mismatch in the output should be a bug in my own compiler.

Csmith has a lot of command line flags to restrict the use of some of the language features. All my tests use:

```c
const char* csmithCmdLine = "csmith.exe "
	"--seed %u "
	"--no-argc "
	"--no-packed-struct "
	"--no-volatiles "
	"--no-volatile-pointers "
	;
```

where `%u` is the seed value in order to be able to reproduce the tests if something fails.

## Debugging

Whenever a test fails with an invalid checksum, I enable `print_hash_value` and diff the outputs. I then put a `printf` for the offending global before and after the call to `func_1()` in `main`. Based on the expected value some bugs were easy to find. If that is not the case, I compile the offending test using `clang-cl` in Visual Studio and put a data breakpoint on the global. This way I narrow the search to the expressions altering the value and try to reproduce it in a minimal example. If that still doesn't work, I've made a version of `safe_math.h` which prints all executed operations and their operands. If after comparing the traces the bug is still not obvious, I put `printf` traces to all `/* block id: xxx */` in the code. Usually after comparing those traces the bug reveals itself. If not, I leave the test and continue with the next one, hoping that it will resolve itself by fixing some other test.

## Bugs found so far

So far I've found several bugs:

- Unnamed 0-width bitfields. 
	- I had an assert for this case but never managed to implement it.
- Parsing of hexadecimal floating point numbers
- Division by constant
	- I had an assert for one specific case which I haven't found a divisor to test it.
	- Codegen for i32 sdiv by negative divisor was wrong.
- Sign-extension and integer promotions of bitfields
- The returned value of a compound assignment expression (`+=`, `-=`, etc) was wrong.
- Bitwise NOT (`~`) integer promotions
- Compound assignment to bitfield member
- Removal of trivial phis with Undefined inputs.
	- Those were removed unconditionally but they can only be removed if the only defined input dominates the phi.
- Don't move `sdiv/srem/udiv/urem` which might trap, outside of loops during LICM.
	- Only ops by constants which are guaranteed to not trap are allowed to be moved.
- Global union initializers
	- The old code expected the initializer to have the same type as the declared type. Unions are special in this regard because the underlying IR type might not match the declared type. Changed the code to allow the initializer to have a different type than the declared type and fill the rest of the allocated space with zeros. Relocations work as expected.

## Other changes made to the compiler

- Added an Undefined value to the compiler.
	- Still a WIP. I also need an `unreachable` instruction to replace branches on undefined conditional values.
- Converted several recursive algorithms (DFS during dominator tree construction, SCCs) to be iterative.
- Changed the AST-to-IR translation code to:
	- have a max recursion depth limit, because Csmith can generate really long expressions.
	- use per-AST node callbacks instead of big switches when generating statements and expressions.
- Generate unwind info during codegen in order to be able to wrap Csmith's main in `__try/__except` blocks and continue testing even if the code crashes.

Still testing...