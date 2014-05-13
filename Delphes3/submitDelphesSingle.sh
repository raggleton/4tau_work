#!/bin/bash
hepmcname="/panfs/panasas01/phys/ra12451/pythia8183/Signal_1prong_500K_bare/signal_1prong_500k_5.hepmc"
rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/Signal_1prong_500K_bare/signal_1prong_500K_5.root"
qsub -v hepmcname=${hepmcname},rootname=${rootname} runDelphesBatch.sh
