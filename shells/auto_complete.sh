#!/bin/bash

_bm_complete() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"

    # Auto-complete para o primeiro argumento: nome do arquivo sem extensão
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        local files=$(ls ../Data/TidyData/ 2>/dev/null | sed -E 's/\.[^.]+$//' | sort -u)
        COMPREPLY=( $(compgen -W "${files}" -- "${cur}") )
        return 0
    fi

    # Auto-complete para o sexto argumento: método
    if [[ ${COMP_CWORD} -eq 6 ]]; then
        local methods="metropolis exact parallel_tempering"
        COMPREPLY=( $(compgen -W "${methods}" -- "${cur}") )
        return 0
    fi
}

# Registra o auto-complete para o script BoltmannMachine.sh
complete -F _bm_complete ./BoltmannMachine.sh