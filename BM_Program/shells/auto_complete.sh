#!/bin/bash

_synt_completions() {
    local cur="${COMP_WORDS[COMP_CWORD]}"
    local opts="exact metropolis parallel_tempering swendsen_wang wang_landau"
    COMPREPLY=( $(compgen -W "${opts}" -- "$cur") )
}

# Registro para o caminho relativo usado para executar o binário
complete -F _synt_completions ../bins/synt
