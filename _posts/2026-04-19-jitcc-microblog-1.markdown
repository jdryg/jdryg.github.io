---
title: "jitcc microblog #1: Function attributes"
date:  2026-04-19 12:00:00 +0200
categories: compilers
tags: attributes
---

## Function Attributes pass

Added an IR pass to calculate function attributes. Currently only `readnone` (aka GCC/Clang `__attribute__((const))`) and `readonly` (aka `__attribute__((pure))`) are supported by the compiler, because these are the only 2 I know how to actually use at the moment. `readnone` helps with GVN and `readonly` helps with MemorySSA. I also have `inline` and `noinline` but those can only be specified in the source code.

I needed this in order to cleanup the IR generated from Csmith tests. All Csmith safe math functions are `readnone` so if their result is unused they can be safely eliminated by DCE.

The code for assigning attributes to a single function is really simple. Start by assuming that the function is `readnone`. Scan the instruction stream from top to bottom and 
- If a `store` is encountered, mark the function as `readwrite` (i.e. no attribute) and stop the scan. 
- If a `load` is encountered, mark the function as `readonly` and continue
- If a `call` is encountered:
  - If it's an indirect `call`, mark the function as `readwrite` and stop the scan.
  - If the called function is marked as `readnone`, ignore the call.
  - If the called function is marked as `readonly`, mark the function as `readonly`.
  - Otherwise, mark the function as `readwrite` and stop the scan.

The issue with calls is that, if the called function hasn't been processed yet, it doesn't have its attributes calculated and the attributes of the current function will be conservative (`readwrite`). In order to avoid this, the pass is a module pass. It's executed on the whole module (with LTO enabled this means it is executed once on the whole program). Having all the functions in a module allows building the call graph, which in turn is used to build Strongly Connected Components. I'm currently using Tarjan's SCC algorithm which returns all SCCs in reverse topological order. Simply by walking the SCC list from head to tail guarantess that a function have its attributes calculated before encountering any `call` to it.

The only issue I haven't solved yet is recursive functions (direct or indirect). Functions in multi-node SCCs (i.e. cycles) are processed in the order they were added to the SCC, which is wrong. It doesn't produce wrong results though. It's just that the assigned attributes are more conservative. As far as I understand, I need an iterative algorithm to calculate attributes of functions participating in a cycle. But then I have to find the entry to the cycle. And if the cycle has multiple entries then things get complicated. I'll leave that for later.
