---
title:  "Generating sets of N variables"
date:   2024-03-24 10:00:00 +0200
categories: code
tags: code
---

Assume we have `N` variables for which we want to generate all sets such that each variable takes `Ki` distinct values. 

For example, if `N=3`, this could be done with 3 nested loops.

```c
double variables[3];
for (uint32_t iVar0 = 0; iVar0 < K[0]; ++iVar0) {
  variables[0] = calcVariableValue(0, iVar0);
  
  for (uint32_t iVar1 = 0; iVar1 < K[1]; ++iVar1) {
    variables[1] = calcVariableValue(1, iVar1);

    for (uint32_t iVar2 = 0; iVar2 < K[2]; ++iVar2) {
      variables[2] = calcVariableValue(2, iVar2);

      // Now do something with this set of values.
      calculate(variables);
    }
  }
}
```

This gets ugly really fast with increasing `N` and it cannot be used when `N` is not known at compile time. Recursion helps in this case.

```c
static void interpolateVariable_r(double* variables, uint32_t N, const uint32_t* K, uint32_t index)
{
  if (index == N) {
    calculate(variables);

    return;
  }

  for (uint32_t iVar = 0; iVar < K[index]; ++iVar) {
    variables[index] = calcVariableValue(index, iVar);
    interpolateVariable_r(variables, N, K, index + 1);
  }
}

double variables[5];
interpolateVariable_r(variables, 5, K, 0);
```

This works fine as long as `N` is small enough to not cause a stack overflow. The only minor issue with this approach is that it's a bit difficult to identify a specific set later. E.g. if the calculation generated an error and you want to find the set with the smallest error.

The following method doesn't use recursion (no stack overflows) and we can easily identify a specific set (and regenerate it later if needed).

```c
uint32_t totalSets = 1;
for (uint32_t index = 0; index < N; ++index) {
  totalSets *= K[index];
}

for (uint32_t iSet = 0; iSet < totalSets; ++iSet) {
  double variables[N];

  uint32_t divisor = 1;
  for (uint32_t i = N; i > 0; --i) {
    const uint32_t index = i - 1;

    const uint32_t iVar = (iSet / divisor) % K[index];
    variables[index] = calcVariableValue(index, iVar);

    divisor *= K[index];
  }

  calculate(variables);
}
```

Let's see an example. Assume `N = 3` and `K[0] = 1`, `K[1] = 5` and `K[2] = 3`. All the sets are shown in the table below.

iSet | iVar[0] | iVar[1] | iVar[2]
-|:-:|:-:|:-:
0 | 0 | 0 | 0
1 | 0 | 0 | 1
2 | 0 | 0 | 2
3 | 0 | 1 | 0
4 | 0 | 1 | 1
5 | 0 | 1 | 2
6 | 0 | 2 | 0
7 | 0 | 2 | 1
8 | 0 | 2 | 2
9 | 0 | 3 | 0
10 | 0 | 3 | 1
11 | 0 | 3 | 2
12 | 0 | 4 | 0
13 | 0 | 4 | 1
14 | 0 | 4 | 2

From the above table it's observed that:
- `iVar[2] = iSet % 3` or `iVar[2] = (iSet / 1) % K[2]`.
- `iVar[1] = (iSet / 3) % 5` or `iVar[1] = (iSet / K[2]) % K[1]`
- `iVar[0] = (iSet / 15) % 1` or `iVar[0] = (iSet / (K[2] * K[1])) % K[0]`


That's all for now. Thanks for reading.
