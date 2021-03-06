# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions
module add shared default-environment
module add tools/git-1.8.4.2
module add languages/python-2.7.6
module add languages/gcc-4.8.1
#module add gcc/4.7.0

export CLICOLOR=1
alias root="root -l"

export PS1="\[\033[36m\][\t]\[\033[m\]:\[\033[32;1m\]\w\[\033[m\] >> "

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/HepMC-2.06.08/x86_64-slc5-gcc43-opt/lib:$HOME/root/root/lib:$HOME/boost_1_55_0_install/lib
export PATH=$PATH:$HOME/root/root/bin:/panfs/panasas01/phys/ra12451/bin

alias duh='du -hcs * | sort -h'
alias l.='ls -d .* --color=auto'
alias la='ls -a'
alias lh='ls -lh'
alias ll='ls -l --color=auto'
alias ls='ls --color=auto'
alias lt='ls -lrth'
alias lS='ls -rSlh'
alias mc='. /usr/libexec/mc/mc-wrapper.sh'
alias pu="ps ux"
alias qstatme='qstat -u ra12451; getNumberInQueue '
alias delEmpty="find . -type f -empty -delete"
alias loadDelphes="gSystem->Load("libDelphes");"
alias Latest="less `ls -rt | tail -n 1`"
alias vi='vim'
#alias vsqueue='qstat | grep veryshort && echo '\''TOTAL: '\'' DAMMIT'
alias vsqueue='qstat | grep veryshort | tee  >(wc -l)'
alias squeue='qstat | grep [^very]short | tee  >(wc -l)'
alias which='alias | /usr/bin/which --tty-only --read-alias --show-dot --show-tilde'


getNumberInQueue(){
	# Put as function, not alias, because bash is a pile of wank and won't parse it correctly without much work
	qstat -u ra12451 | awk 'BEGIN{n=0; nR=0; nQ=0; nC=0;}{if($2=="ra12451") n++; if($10=="R") nR++; else if($10=="Q") nQ++; else if($10=="C") nC++;}END{print "Jobs queueing: ",nQ; print "Jobs running: ",nR; print "Jobs completed: ", nC; print "Total jobs: ",n;}'
}
