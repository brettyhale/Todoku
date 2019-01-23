## Todoku
__Todoku__ is a <i>general</i> Sudoku solver for <b>NxN x NxN</b> puzzles (with an eternal 'TODO' list...)

There are limits placed on <b>N</b> by the implementation. The current range for <b>N</b> is : <b>1</b> to <b>6</b></br>
Character encodings are as follows:

`N = 1` : '<b>0</b>' or '<b>.</b>' denote an unsolved cell. The trivial solution is: '<b>1</b>'</br>

`N = 2` : '<b>0</b>' or '<b>.</b>' denote an unsolved cell. Solution values are: '<b>1</b>' .. '<b>4</b>' - (4 x 4 Primer for kids) (16 cells)</br>

`N = 3` : '<b>0</b>' or '<b>.</b>' denote an unsolved cell. Solution values are: '<b>1</b>' .. '<b>9</b>' - (Standard Sudoku) (81 cells)</br>

`N = 4` : '<b>.</b>' denotes an unsolved cell. Solution values are: '<b>0</b> .. '<b>9</b>', '<b>A</b>' .. '<b>F</b>' - (16 x 16 Hexadoku) (256 cells)</br>

`N = 5` : '<b>0</b>' or '<b>.</b>' denote an unsolved cell. Solution values are: '<b>A</b>' .. '<b>Y</b>' - (25 x 25 Alphadoku) (625 cells)</br>

`N = 6` : '<b>.</b>' denotes an unsolved cell. Solution values are: '<b>0</b> .. '<b>9</b>', '<b>A</b>' .. '<b>Z</b>' - (36 x 36 Alphanumeric) (1296 cells)</br>

---

#### Implementation notes:

Todoku is distributed as a single C++ source file: `todoku.cc`. It has been an exercise in using various C++11 language features (requiring a C++11 compiler), concise formatting, and pertinent commentary. It is also an experiment in trying to achieve a balance between recursive algorithms and heuristic strategies.

The test vectors (datasets) are not mine. I've preserved commentary in these files where practical, but cannot - in general - vouch for their public domain status. I'm happy to amend, replace, (or even remove) datasets at the request of their original owners. On the subject of test vectors, I feel that much of the <b>todoku</b> project code could be adapted to <i>generate</i> puzzles. This would require considerable re-factoring, as well as making the `Game` class somewhat more robust to prevent illegal states, etc. I have not investigated efficient methods of puzzle generation that guarantee unique solutions.

Hard limits for '<b>N</b>' outside the I/O functionality (i.e., character <-> value mapping) can obviously be much higher. The `unsigned int` type is used extensively. This has a mandated minimum of 16 bits; but is a 32 bit type on almost all relevant ABIs for general purpose computing at this time. I have not performed anything approaching a formal verification on a maximum '<b>N</b>' size, but I'm confident that `N = 15` (a staggering `50625` cells!) can be implemented safely - ignoring the I/O code.

Of course, Sudoku is a known <b>NP-complete</b> problem, so <b>N</b> can always be increased to defeat any amount of (classical) computing power... regardless of the various algorithms and strategies applied.
