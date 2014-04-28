# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# User specific aliases and functions
module add shared default-environment
module add tools/git-1.8.4.2
module add languages/python-2.7.6
module add gcc/4.7.0

export CLICOLOR=1
alias root="root -l"

export PS1="\[\033[36m\][\t]\[\033[m\]:\[\033[32;1m\]\w\[\033[m\] >> "

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/panfs/panasas01/phys/ra12451/HepMC-2.06.08/x86_64-slc5-gcc43-opt/lib:/panfs/panasas01/phys/ra12451/root/root/lib:/panfs/panasas01/phys/ra12451/boost_1_55_0_install/lib
export PATH=$PATH:/panfs/panasas01/phys/ra12451/root/root/bin:/panfs/panasas01/phys/ra12451/bin

alias duh='du -hcs * | sort -h'
alias l.='ls -d .* --color=auto'
alias la='ls -a'
alias lh='ls -lh'
alias ll='ls -l --color=auto'
alias ls='ls --color=auto'
alias lt='ls -lrth'
alias mc='. /usr/libexec/mc/mc-wrapper.sh'
alias qstatme='qstat -u ra12451'
alias vi='vim'
alias vsqueue='qstat | grep veryshort && echo '\''TOTAL: '\'' 545'
alias which='alias | /usr/bin/which --tty-only --read-alias --show-dot --show-tilde'
