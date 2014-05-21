#!/bin/bash

for i in {31..50}
do
	echo $i
	./DelphesHepMC examples/delphes_card_CMS_bare.tcl QCDc_mu_pthatmin20_Mu17_Mu8_bare/QCDc_mu_pthatmin20_Mu17_Mu8_${i}.root ~/pythia8183/QCDc_mu_pthatmin20_Mu17_Mu8/qcdc_pthatmin20_${i}_HLT.hepmc
done
