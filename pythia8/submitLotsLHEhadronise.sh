for file in ${HOME}/pythia8185/Signal_1prong_500K_bare/GG_H_aa_8_4taus_decay_500K_[0-9]*.lhe
do
	echo "Doing $file"
	filename=${file##*/}
	out=${filename%.lhe}
	number=${out%-single}
	echo $number
	number=${number#GG_H_aa_8_4taus_decay_500K_}
	qsub -N ma8_${number} -v LHEname=${file},hepmcname="${HOME}/pythia8185/Signal_1prong_500K_bare/"${out%-single} runLHEhadronise.sh
	sleep 10
done
