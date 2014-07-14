#!/bin/bash

# Script to do lots of Delphes submissions locally

# for name in ${HOME}/pythia8183/examples/qcdb_*.hepmc
# for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20/qcdb_pthatmin20_*.hepmc
# for name in ${HOME}/pythia8183/QCDb_mu_pthatmin20_Mu17_Mu8/qcdb_pthatmin20_Mu17_Mu8_*.hepmc
# for name in ${HOME}/pythia8183/QCDc_mu_pthatmin20_Mu17_Mu8/qcdc_pthatmin20_*.hepmc
# for name in ${HOME}/pythia8183/QCDScatter_mu_pthatmin20_Mu17_Mu8/qcdScatter_pthatmin20_*.hepmc
# for name in ${HOME}/pythia8183/QCDAll_mu_pthatmin20_Mu17_Mu8/qcdAll_pthatmin20_*.hepmc
# for name in ${HOME}/pythia8183/QCDAll_NEW_mu_pthatmin20_Mu17_Mu8/qcdAll_pthatmin20_NEW_*.hepmc
# for name in ${HOME}/pythia8185/QCDb_HLT/b_5000_*.hepmc
# for name in ${HOME}/pythia8185/QCDb_HLT/b_10000_*.hepmc
# hepmcname="scatterbc_250_"
hepmcname="GG_H_aa_8_4taus_decay_500K_"
# for name in ${HOME}/pythia8185/QCDScatter_HLT/${hepmcname}*.hepmc
for name in ${HOME}/pythia8183/Signal_1prong_500K_bare/${hepmcname}*_HLT.hepmc
do
	filename=${name##*/}
	# num=${filename##*_}
	# pureNum=${num%.hepmc}
	# stem=${filename%_NoHLT.hepmc}
	# stem=${filename%_HLT.hepmc}
	stem=${filename%-single_HLT.hepmc}
	# pureNum=${stem#qcdb_pthatmin20_Mu17_Mu8_}
	# pureNum=${stem#qcdc_pthatmin20_}
	# pureNum=${stem#qcdScatter_pthatmin20_}
	# pureNum=${stem#qcdAll_pthatmin20_}
	# pureNum=${stem#qcdAll_pthatmin20_NEW_}
	# pureNum=${stem#b_5000_}
	# pureNum=${stem#b_10000_}
	pureNum=${stem#$hepmcname}
	#if [ "$pureNum" -gt "220" ]
	#then
		echo "Doing "$filename
		# out="QCDc_mu_pthatmin20_Mu17_Mu8_$pureNum"	
		# out="QCDScatter_mu_pthatmin20_Mu17_Mu8_$pureNum"	
		# out="QCDAll_mu_pthatmin20_Mu17_Mu8_$pureNum"	
		# out="QCDAll_NEW_mu_pthatmin20_Mu17_Mu8_$pureNum"	
		# out="QCDb_HLT_$pureNum"	
		# out="QCDbcScatter_HLT_250_$pureNum"	
		out="Signal_HLT_$pureNum"	
		# folder="QCDScatter_mu_pthatmin20_Mu17_Mu8_bare"
		# folder="QCDAll_mu_pthatmin20_Mu17_Mu8_bare"
		# folder="QCDAll_NEW_mu_pthatmin20_Mu17_Mu8_bare"
		# folder="QCDbcScatter_HLT_bare"
		folder="Signal_1prong_HLT_bare"
		echo "Writing to "${folder}"/"${out}.root
		# echo $out
		# echo $name
		# qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDb_mu_pthatmin20_Mu17_Mu8_bare/"${out}.root runDelphesBatch.sh
		# qsub -v hepmcname=${name},rootname="/panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDc_mu_pthatmin20_Mu17_Mu8_bare/"${out}.root runDelphesBatch.sh
		# sleep 10
		# ../DelphesHepMC ../examples/delphes_card_CMS_bare.tcl /panfs/panasas01/phys/ra12451/Delphes-3.0.12/QCDc_mu_pthatmin20_Mu17_Mu8_bare/${out}.root $name
		../DelphesHepMC ../examples/delphes_card_CMS_bare.tcl /panfs/panasas01/phys/ra12451/Delphes-3.1.2/${folder}/${out}.root $name
	#fi
done
