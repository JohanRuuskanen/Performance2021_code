### Reproducing the examples/evaluation

The code in this repository can be used to reproduce the results seen in the following paper available [here](https://doi.org/10.1016/j.peva.2021.102231)

> Johan Ruuskanen, Tommi Berner, Karl-Erik Årzén, Anton Cervin, *Improving the mean-field fluid model of processor sharing queueing networks for dynamic performance models in cloud computing*. Performance Evaluation, Volume 151, 2021

The implementation is done in `Julia`, but requires the two following external software tools to run

​	`LINE`, http://line-solver.sourceforge.net/
​	`JMT`, http://jmt.sourceforge.net/

which are supplied in the repository. The requirements are listed below, with the versions used during development and testing

   - `Julia v1.6.1`
   - `MATLAB R2018a` (for LINE)
   - `Java openjdk 11.0.11` (for JMT)
   - `Python 3.6.9` with `matplotlib 3.1.1` (for plotting)

To run, start the `Julia` REPL in this directory and type

```julia
] activate .
] instantiate
] precompile
```

to install the required julia packages. The files for running the evaluation can then be found in the `manuscript/evaluation` folder, while the examples can be found in any of the `manuscript/example*` folders. They should be runnable out-of-the-box, after changing the `basepath` variable. 









