---
title: "Local Value Numbering"
date:  2026-01-07 19:45:00 +0200
categories: compilers
tags: optimizations, lvn
---

I've been writing a JIT C compiler for the past 8 months and I thought it might be good to start documenting some of my experiments. I'm in no way suitable for talking about compilers or optimizations, so take the following with a large grain of salt. 

### Intermediate Representation

The IR I'll be using in this post is a 3-address code (TAC) IR. I call it LIR in the compiler (Low-level IR) but it can be considered a middle-level IR. It includes high-level constructs to hide ABI-specific details (va_start/va_end/etc instructions and specific function argument operands), but it also supports complex addressing modes like x64 (i.e. `[%vr0 + %vr1 * 8 + 16]`). There are also a High-level IR (similar to LLVM) and a Machine-level IR (identical to x64 assembly with both hardware and virtual registers).

There are 8 kinds of operands in LIR:
- Register 
    - Only virtual registers. There is no concept of hardware registers at this level
    - There are 2 classes of registers: Integer and Floating point, each with their own types.
- Constant (either integer or floating point)
- Basic Block 
    - Used only by branch and phi instructions
- Generic Memory Reference (`[base_reg + index_reg * scale + displacement]`)
    - Both `base_reg` and `index_reg` can be set to `None`
    - Both `base_reg` and `index_reg` are expected to be integer registers
- External Symbol 
    - RIP-relative memory references used to point to global variables or functions.
- Stack Object
    - Stack-relative memory references used for function-local variables.
- Comparison Predicate
    - Used only by comparison and conditional branch instructions.
- Function Argument
    - Used only by the `copy` instruction to move function arguments to a virtual register.

It's a load/store IR, meaning that most instructions accept only registers as operands, some accept constants and only specific instructions accept memory operands. 

`call` is a multi-operand instruction, breaking the strict 3-address code concept. It's the only instruction with an optional (i.e. `NULL`) destination operand (all other instructions either always have a valid destination operand or don't have any destination operand at all). All call arguments are always registers.

There is a `phi` instruction to support code in SSA form, which is mandatory for some transformation passes, but the IR API does not impose any restriction on that. In other words, a transformation pass is free to define a virtual register multiple times if it needs to. `phi` is also a multi-operand instruction, with a single destination operand and 1 pair of (register, basic block) operands for each of the basic block's predecessors.

### High-level algorithm

Assuming the code is in SSA form, the basic high-level algorithm for hash-based Local Value Numbering (LVN) is:

```
for each basic block BB {
    clear hash table T
    clear value number to register map

    for each instruction INSTR in BB {
        if INSTR is in the form DST = opcode SRC1, SRC2, ... {
            HASH = calculate instruction hash
            if HASH is in table T {
                (VN, PREV) = get value number and instruction of HASH from T
                map VN to DST
                map DST to the destination reg of PREV
                mark INSTR as dead
            } else {
                VN = generate new value number
                insert (VN, INSTR) to T
                map VN to DST
                map DST to DST
            }
        }
    }
}
```

The output of the above algorithm is a mapping between virtual registers for the whole function. This mapping is used to rewrite instructions in a second pass. Finally all dead instructions are removed from the function.

### Calls

As mentioned above `call` is a special instruction from the IR's perspective but it's also a special instruction from the LVN's perspective. Even though its form is `DST = opcode SRC1, SRC2, ...` as expected by the algorithm above, `DST` might be `NULL` but more importantly we don't want to hash and eliminate calls, because this will change the behavior of the program. 

E.g. multiple `%vr_i = call myFunc %vr0` should not be replaced by the first `call`'s result because the algorithm does not know if `myFunc` has side-effects or touches global state. If the algorithm knew more about `myFunc` (e.g. leaf function without side effects) it might be possible to value number those calls and eliminate all but the first, but since that info is not (currently) available, `call` and call-like instructions are handled separately in the inner loop.

```
    for each instruction INSTR in BB {
        if isCall(INSTR) {
            if destination operand DST is not NULL {
                VN = generate new value number
                map VN to DST
                map DST to DST
            }
        } else if INSTR is in the form DST = opcode SRC1, SRC2, ... {
            // Same as before...
        }
    }
```

`isCall` checks if the instruction's opcode is `call`-like. Currently it checks for `call`, `va_start`, `va_copy` and `va_end`. All those are treated as function calls.

### Phis

The other special instruction is `phi`. `phi`s always have a destination operand (register). That register always get a new value number. I haven't tried hashing phis because I don't expect to end up with multiple identical phis in the same basic block. Also I don't know if it's valid to eliminate phis because this might complicate the rest of the code.

```
    for each instruction INSTR in BB {
        if isCall(INSTR) {
            // Same as before...
        } else if INSTR is DST = phi [val1, BB1], [val2, BB2], ... {
            VN = generate new value number
            map VN to DST
            map DST to DST
        } else if INSTR is in the form DST = opcode SRC1, SRC2, ... {
            // Same as before...
        }
    }
```

### Terminators

Terminator instructions (branchs and `ret`) are ignored by the algorithm. They are by definition not in the expected form to be hashed and there is no point to, somehow, value-numbering them.

### Loads and stores

Instructions which load data from memory should be value numbered because theoretically there is no benefit in loading the same value from memory multiple times if it hasn't changed. They are also in the correct form to be hashed. The issue with loads is stores. Hashing a `load` instruction based only on its operands will give wrong results if there is an intervening store between 2 identical loads.

```
%vr1 = load dword [%vr0 + 8];
store dword [%vr2], %vr1;
%vr3 = load dword [%vr0 + 8];
; Other arithmetic instructions and loads
%vr4 = load dword [%vr0 + 8];
```

Without alias analysis the algorithm cannot know if writing 4 bytes to `[%vr2]` overwrites the value loaded from `[%vr0 + 8]`. The algorithm needs a way to identify that the 2nd `load` is not redundant and it should assign `%vr3` a different value number than `%vr1`. But it also needs to figure out that the 3rd load is redundant because there is no intervening instruction which affects the already loaded value into `%vr3`. 

The easiest way to handle this is to introduce a new global value number for memory (`mem_vn`). Whenever the algorithm encounters an instruction which might affect the memory (`store`, `call`, etc), `mem_vn` is assigned a new value number. When hashing `load` instructions `mem_vn` is included in the hashed data and as a result a different hash is produced for each load, from the same memory location, if `mem_vn` changed at some point.

### Memory objects

As mentioned in the beginning there are 3 kinds of memory-related operands. Generic memory references, which follow x64 addressing modes, stack objects and external symbols. Stack objects are disjoint, even though they actually end up sharing the same stack. The same is true for external symbols, even though they most probably end up sharing the same section in the final image. 

Assuming the initial program code was correct (e.g. all accesses to a specific stack object are within its allocated space), we can take the value numbering of stack objects and external symbols one step further by assigning different value numbers to different segments based on how they are accessed. E.g.

```
%vr1 = load dword [stack_obj.0];
store dword [stack_obj.1 + 4], %vr1;
%vr2 = load dword [stack_obj.0];
%vr3 = load dword [stack_obj.1 + 4];
```

In this case, the 2nd load is redundant because the `store` to `[stack_obj.1 + 4]` cannot overwrite the value loaded into `%vr1` from `[stack_obj.0]`. In other words, the `store` does not affect the value number assigned to the [0, 4) segment of `stack_obj.0`. Additionally, it might help if the 3rd load was eliminated by forwarding `%vr1` (i.e. the last stored value number to `[stack_obj.1 + 4]`) to `%vr3`.

This is implemented by creating a memory object for each stack object and external symbol. This memory object is split into segments, based on how it's accessed, and each segment keeps its own "memory VN" and its "last stored VN". Initially all memory objects get a single continuous segment [0, UINT32_MAX], its memory VN gets a new value and last stored VN is set to an invalid value. Even though the code knows the true size of stack objects are this point, this might not be true for external symbols defined in separate translation units.

```
typedef struct mem_segment_t
{
    mem_segment_t* m_Next;
    mem_segment_t* m_Prev;
    uint32_t m_Offset;       // Segment's offset
    uint32_t m_Size;         // Segment's size
    uint32_t m_MemVN;        // mem_vn of the segment used when hashing loads
    uint32_t m_LastStoredVN; // VN of the last stored reg to this segment; might be invalid
} mem_segment_t;

typedef struct mem_object_t
{
    mem_segment_t* m_SegmentListHead;
    mem_segment_t* m_SegmentListTail;
} mem_object_t;
```

All memory objects are kept into a hashmap. Every time a `store` to a generic memory reference is encountered, the hashmap is cleared and all information about memory objects is discarded. Store-to-load forwarding is implemented in case of `load`s just before hashing the instruction. If the forwarding succeeds, hashing is skipped and the algorithm moves on to the next instruction. Otherwise, the `load` is hashed as a regular instruction, including all the VNs of the segments covering the accessed part.

Let's see how each instruction from the above snippet is processed.
- `%vr1 = load dword [stack_obj.0]`
    - This is a `load`. Before hashing, check if there is a valid VN stored to `dword [stack_obj.0]`.
    - Stack object #0 has no segment at offset 0 with size 4 yet. It is split into 2 segments: [0, 4) + [4, UINT32_MAX].
    - Both new segments get the same VNs as the parent segment (i.e. invalid last stored VN).
    - Since the new segment's last stored VN is invalid, there is nothing to do at this point.
    - The `load` is hashed using segment's [0, 4) VN.
    - Nothing is found in the hashtable for this hash. A new VN is assigned to `%vr1` and the `load` is inserted into the table.
- `store dword [stack_obj.1 + 4], %vr1`
    - This is a `store`. The source register's VN (`%vr1`) should be stored to the corresponding memory segment.
    - Stack object #1 is split into 3 segments: [0, 4) + [4, 8) + [8, UINT32_MAX]. 
    - All three get the same value numbers as the initial segment.
    - Segment's [4, 8) last stored VN is set to `%vr1`'s VN and its VN gets a new value.
    - Stack object #0 VNs are not affected.
    - No hashing is performed because `store` does not have a destination register.
- `%vr2 = load dword [stack_obj.0]`
    - This is a `load`. Before hashing, check if there is a valid VN stored to `dword [stack_obj.0]`.
    - A segment [0, 4) for stack object #0 already exists so no splitting is necessary. 
    - Its last stored VN is invalid since we haven't encountered a `store` to this segment yet.
    - The `load` is hashed using the segment's [0, 4) VN. Since the VN hasn't changed, the hash ends up being the same as the one calculated for the first instruction and as a result `%vr2` is mapped to `%vr1` and the `load` is marked as dead.
- `%vr3 = load dword [stack_obj.1 + 4]`
    - This is a `load`. Before hashing, check if there is a valid VN stored to `dword [stack_obj.1 + 4]`.
    - A segment [4, 8) for stack object #1 already exists. Its last stored VN is equal to `%vr1`'s VN.
    - As a result `%vr3` is mapped to `%vr1` and the `load` is marked as dead.
    - No hashing is performed.

Stores and loads to/from external symbols are handled in the same way. Doing store-to-load forwarding this way for generic memory references works only if the global `mem_vn` hasn't changed and only for the last generic memory store (i.e. keep track of the base/index/scale/displacement of the last generic store and if a load matches those, use the stored VN).

### Commutative operations

Commutative operations (`add`, `mul`, `fadd`, etc) should be hashed using some predefined ordering of their operands. Currently, all arithmetic instructions support either (register, register) or (register, constant) operands. When both operands are registers they are first ordered based on their VNs and hashed in that order. This allows calculating the same hash for both `add %vr1, %vr2` and `add %vr2, %vr1`.

### Comparisons

Comparisons can be handled in the same way as the commutative operations. The only change required is to invert the comparison predicate operand. E.g. a signed less-equal comparison `icmp sle, %vr2, %vr1` might be hashed as a signed greater-equal `icmp sge, %vr1, %vr2` due to the order of register VNs. 

### The end

That's it. If I haven't forgot any major detail, that's how LVN is currently implemented in my SSA TAC IR. Thanks for reading.
