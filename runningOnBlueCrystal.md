# Running PYTHIA & Delphes on batch system

(e.g Blue Crystal)

## THINGS TO INSTALL 

(in this order)

- ROOT
- HEPMC
- PYTHIA8
- Delphes 3



To link to HEPMC:

- Download pre built hepmc
- Compile pythia with hepmc options (see README.hepmc)
- When submitting, make sure you add 

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/panfs/panasas01/phys/ra12451/HepMC-2.06.08/x86_64-slc5-gcc43-opt/lib
```

to `~/.bashrc`, otherwise will complain about the libhepmc 

**How does qsub know where libraries etc are?** -> it uses your `PATH` and `LD_LIBRARY_PATH` s to import the correct libs etc

## How to use variable in shell script passed to `qsub`

```
qsub -v param_name=value script_name
```

in your `script_name`

e.g.

```
echo ${param_name}
```

see: https://stackoverflow.com/questions/3504081/parameter-for-shell-scripts-that-is-started-with-qsub
