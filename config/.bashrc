# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
if [ "$TERM" = "xterm" ] || [ "$TERM" = "xterm-256color" ]
then
    #we're on the system console or maybe telnetting in
  # export PS1="\[\033]2;\h:\u \w\007\033[33;1m\]\u \[\033[36;1m\]\w\[\033[0m\]\n\[\e[32;1m\]$ \[\e[0m\]"
    #export PS1="\[\033]2;\h:\u \$PWD\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]\$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
    export PS1="\[\033]2;\h:\u \$PWD\007\033[33;1m\]\u@\h \033[35;1m\t \e[4;1m \$(date +%d.%m.%Y) \n\033[0m\[\033[36;1m\]\$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
    #export PS1="\[\033]2;\h:\u \w\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]\w\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
else
    #we're not on the console, assume an xterm
    #export PS1="\[\033]2;\h:\u \$PWD\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]\$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
    export PS1="\[\033]2;\h:\u \$PWD\007\033[33;1m\]\u@\h \033[35;1m\t \e[4;1m \$(date +%d.%m.%Y) \n\033[0m\[\033[36;1m\]\$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
    #export PS1="\[\033[33;1m\]\h:\u \[\033[37;1m\]\w\$\[\033[0m\]"
fi
export CLICOLOR=1
export LSCOLORS=ExFxBxDxCxegedabagacad

export HISTSIZE=200000000
export HISTFILESIZE=500000000


export PATH=/ds3200_1/users_root/yitingshuang/applications/bowtie2-2.3.4-linux-x86_64:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/SPAdes-3.12.0-linux/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/OGA/scripts:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/PGA:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/corset-1.07-linux64:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/Salmon-0.8.2_linux_x86_64/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/transrate-1.0.4/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/TransDecoder-TransDecoder-v5.3.0:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/TransDecoder-TransDecoder-v5.3.0/util:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/Digest-MD5-2.52/usr/lib64/perl5/vendor_perl/Digest:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/plot-ks-master:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/KaKs_Calculator2.0/bin/Linux:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/TreeShrink/:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/OrthoFinder/fastme-2.1.5/binaries/:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/OrthoFinder/mmseqs2/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/OrthoFinder/Diamond/:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/OrthoFinder/OrthoFinder-2.3.3_source/orthofinder/:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/applications/MAPS/barkerlab-maps/:$PATH
export PATH=//ds3200_1/users_root/yitingshuang/applications/ASTRAL-master/:$PATH













#PATH="/ds3200_1/users_root/yitingshuang/perl5/bin${PATH:+:${PATH}}"; export PATH;
#PERL5LIB="/ds3200_1/users_root/yitingshuang/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;
#PERL_LOCAL_LIB_ROOT="/ds3200_1/users_root/yitingshuang/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"; export PERL_LOCAL_LIB_ROOT;
#PERL_MB_OPT="--install_base \"/ds3200_1/users_root/yitingshuang/perl5\""; export PERL_MB_OPT;
#PERL_MM_OPT="INSTALL_BASE=/ds3200_1/users_root/yitingshuang/perl5"; export PERL_MM_OPT;

#####################################
####LH add start
alias tma='~/lh/anaconda3.1/bin/tmux a'
alias bj='bjobs'
alias bl='bjobs |less'
alias bgrep='bjobs | grep'
alias bjt='bjobs|tail'
alias bjt='bjobs|les'
alias ls='ls --color=tty'
alias bs='bash'
alias les='less -SN'
alias lese='less -SN'
alias ll='ls -hl'                              # long list
alias lt='ls -lth'
alias lth='ls -lth |head'
#alias rm='rm -i'
alias vi='vim'
alias vc='vim ~/.bashrc'
alias sc='source ~/.bashrc'
alias sfat='ssh -Y fat01'
alias s1='ssh -Y cn01'
alias s2='ssh -Y cn02'
alias s3='ssh -Y cn03'
alias s4='ssh -Y cn04'
alias chd='chmod 755 *.sh *.pl *.py'
alias qt='qstat -u yitingshuang|grep "yitingshuang"'
alias tp='top -u yitingshuang'
alias addpwd='echo "export PATH=$PWD:\$PATH" >>~/.bashrc'
alias toadmin='ssh admin@login01'
alias tog='ssh galaxy_user@login01'
alias to1='ssh cn01'
alias to2='ssh cn02'
alias to3='ssh cn03'
alias to4='ssh cn04'
alias big='ssh fat01'
#alias qs='qsub  -V -b y -N output -cwd -l h_vmem=1G'
alias qs='qsub  -V -b y -N output -cwd '
alias ql='qstat -u "*" -f'
alias bsub512='bsub -q Q104C512G_X4  -o output.%J -e error.%J'
alias bsub1T='bsub -q Q64C1T_X4  -o output.%J -e error.%J'
alias bsub2T='bsub -q Q48C2T_X1  -o output.%J -e error.%J'
alias bsub3T='bsub -q Q64C3T_X1  -o output.%J -e error.%J'
alias bsubI='bsub -q Q104C512G_X4  -Ip bash'
alias bsubI1T='bsub -q Q64C1T_X4 -Ip bash'

alias newbsub='echo  "#!/bin/bash
#BSUB -J hybrid_job_name      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail 

ROOT=$PWD"'

body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}

tt () {
    perl /ds3200_1/users_root/yitingshuang/lh/bin/lh_bin/print_head.pl $1 |column -s $'\t' -t   |less -S
}

t () {
    column -s $'\t' -t  "$1" | less -S
}
tailf() ( # args: <file> [<number-of-header-lines>]
  trap 'tput csr 0 "$((LINES-1))"' INT
  tput csr "$((1+${2-1}))" "$((LINES-1))"
  tput clear
  {
    head -n"${2-1}"
    printf "%${COLUMNS}s\n" "" | tr ' ' =
    tail -n "$((LINES-1-${2-1}))" -f
  } < "$1"
)



# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ds3200_1/users_root/yitingshuang/lh/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ds3200_1/users_root/yitingshuang/lh/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/ds3200_1/users_root/yitingshuang/lh/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/ds3200_1/users_root/yitingshuang/lh/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
# export PATH="/ds3200_1/users_root/yitingshuang/lh/anaconda3/bin:$PATH"  # commented out by conda initialize
export HB="/ds3200_1/users_root/yitingshuang/lh/bin/bundle/lh_bin"
export BD="/ds3200_1/users_root/yitingshuang/lh/bin/bundle"
alias vw="vim ${BD}/tools_wrapper.md"

export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/sratoolkit.2.9.0-centos_linux64/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/bundle/lh_bin:$PATH

#ActivePerl
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/ActivePerl-5.24/bin:$PATH
export PERL5LIB=/ds3200_1/users_root/yitingshuang/lh/bin/ActivePerl-5.24/lib:$PERL5LIB
export MANPATH=/ds3200_1/users_root/yitingshuang/lh/bin/ActivePerl-5.24/man:$MANPATH


# added by Anaconda3 installer
# . /ds3200_1/users_root/yitingshuang/lh/anaconda3/etc/profile.d/conda.sh  # commented out by conda initialize
#$. /ds3200_1/users_root/yitingshuang/lh/anaconda2/etc/profile.d/conda.sh
#$export PATH=/ds3200_1/users_root/yitingshuang/lh/anaconda2/bin:$PATH

export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/gatk-4.1.4.0:$PATH

#Autojump
[[ -s /ds3200_1/users_root/yitingshuang/.autojump/etc/profile.d/autojump.sh ]] && source /ds3200_1/users_root/yitingshuang/.autojump/etc/profile.d/autojump.sh

export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/rar:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/SF2:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/.aspera/connect/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/bundle:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/XPCLR/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/selscan:$PATH

#VCFtools
export PERL5LIB=/ds3200_1/users_root/yitingshuang/lh/bin/vcftools/lib:$PERL5LIB
export MANPATH=/ds3200_1/users_root/yitingshuang/lh/bin/vcftools/share/man:$MANPATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/vcftools/bin:$PATH

#boost package
export LDFLAGS=-L/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0/lib
export CFLAGS=-I/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0/include
export CXXFLAGS=-I/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0/include
export CPPFLAGS=-I/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0/include
export LD_LIBRARY_PATH="/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0/lib":$LD_LIBRARY_PATH
export BOOST_ROOT=/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0
export BOOST_LIBRARYDIR=/ds3200_1/users_root/yitingshuang/lh/lib/boost_1_66_0/lib
#export CFLAGS=-I/ds3200_1/users_root/yitingshuang/lh/anaconda3/include/python3.6m/
#export CXXFLAGS=-I/ds3200_1/users_root/yitingshuang/lh/anaconda3/include/python3.6m/
#export CPPFLAGS=-I/ds3200_1/users_root/yitingshuang/lh/anaconda3/include/python3.6m/
export LD_LIBRARY_PATH=/ds3200_1/users_root/yitingshuang/lh/lib/lib64:$LD_LIBRARY_PATH

#bamtools
export LDFLAGS="-L/ds3200_1/users_root/yitingshuang/lh/lib/lib64/ "$LD_FLAGS
export CFLAGS="-I/ds3200_1/users_root/yitingshuang/lh/lib/include/bamtools "$CFLAGS
export CXXFLAGS="-I/ds3200_1/users_root/yitingshuang/lh/lib/include/bamtools "$CXXFLAGS
export CPPFLAGS="-I/ds3200_1/users_root/yitingshuang/lh/lib/include/bamtools "$CPPFLAGS

#common lib
export CFLAGS="-I/ds3200_1/users_root/yitingshuang/lh/lib/include "$CFLAGS
export CXXFLAGS="-I/ds3200_1/users_root/yitingshuang/lh/lib/include "$CXXFLAGS
export CPPFLAGS="-I/ds3200_1/users_root/yitingshuang/lh/lib/include "$CPPFLAGS

export SE_HOME=/ds3200_1/users_root/yitingshuang/lh/bin/snpEff
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/wtdbg2/wtdbg2:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/wtdbg2/minimap2-2.17_x64-linux:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/RepeatMasker:$PATH
export PERL5LIB=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/RepeatMasker:$PERL5LIB
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker/bin:$PATH
#Below is Maker v3.01.02
#export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker/bin:$PATH

#augustus
#export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus/bin:$PATH
#export AUGUSTUS_CONFIG_PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus/config
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config

#snap
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap:$PATH
export ZOE=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe

export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/tRNAscan-SE-1.3.1:$PATH
export PERL5LIB=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/tRNAscan-SE-1.3.1/tRNAscanSE:$PERL5LIB
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nanopolish/nanopolish:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/RepeatModeler/ncbi-rmblastn-2.2.28/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/RepeatModeler/RECON-1.08/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/RepeatModeler/nseg:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/RepeatModeler/RepeatScout-1:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/RepeatModeler/RepeatModeler-open-1.0.11:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/NextDenovo:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/NextDenovo/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/fastx_toolkit:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/FastUniq/source:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/RepeatMasker/util:$PATH

# The next line updates PATH for the Google Cloud SDK.
if [ -f '/ds3200_1/users_root/yitingshuang/lh/bin/google-cloud-sdk-239.0.0/path.bash.inc' ]; then . '/ds3200_1/users_root/yitingshuang/lh/bin/google-cloud-sdk-239.0.0/path.bash.inc'; fi

# The next line enables shell command completion for gcloud.
if [ -f '/ds3200_1/users_root/yitingshuang/lh/bin/google-cloud-sdk-239.0.0/completion.bash.inc' ]; then . '/ds3200_1/users_root/yitingshuang/lh/bin/google-cloud-sdk-239.0.0/completion.bash.inc'; fi
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/last/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/texlive/x86_64-linux:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/p7zip/bin:$PATH

##Latex
export MANPATH=/ds3200_1/users_root/yitingshuang/lh/bin/texlive/2018/texmf-dist/doc/man:$MANPATH
export INFOPATH=/ds3200_1/users_root/yitingshuang/lh/bin/texlive/2018/texmf-dist/doc/info:$INFOPATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/texlive/2018/bin/x86_64-linux:$PATH

##Genblast
#Use this env_var to find formatdb and genblast
export GBLAST_PATH=/ds3200_1/users_root/yitingshuang/lh/bin/genblast_v139

####LH add end
#####################################
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/htop:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/allHiC/bin/ALLHiC/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/allHiC/bin/ALLHiC/scripts:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/lastz-1.04.03/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/bin/genblast_v139:$PATH