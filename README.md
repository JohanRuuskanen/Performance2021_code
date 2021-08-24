### Instructions for running the evaluation/example code

The supplied code includes the following two external software tools

​	`LINE`, http://line-solver.sourceforge.net/
​	`JMT`, http://jmt.sourceforge.net/

Requires `Julia`, `MATLAB` (for `LINE`) and `Java` (for `JMT`) to run, and `Python` with `matplotlib` for plotting. 

Tested with:
   - `Julia v1.5.0`
   - `MATLAB R2018a`
   - `Java openjdk 11.0.11`
   - `Python 3.6.9`
           - `matplotlib 3.1.1`

To run, start the `Julia` REPL in this directory and type

```julia
] activate .
] instantiate
```

to install the required dependencies. The files for running the evaluation can then be found in the `manuscript/evaluation` folder, while the examples can be found in any of the `manuscript/example*` folders. They should be runnable out-of-the-box, after changing the `basepath` variable. 