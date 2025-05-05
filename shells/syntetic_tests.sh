#!/bin/bash

# -------------------------------
# CONFIGURAÇÕES INICIAIS
# -------------------------------

# Lista de tamanhos de spins
N_SPINS=(20 40 60 80)

# Número de execuções por N_spins
N_RUNS=20

# Diretório base para dados e logs (caminho absoluto)
DATA_DIR="$(pwd)/../data_syntetic"
mkdir -p "$DATA_DIR"

ARGS_FILE="$DATA_DIR/input_gen_data.txt"
rm -f "$ARGS_FILE"

LOG_DIR="$DATA_DIR/logs_gen_data"
mkdir -p "$LOG_DIR"

# -------------------------------
# CAMINHO ABSOLUTO PARA O GEN_DATA
# -------------------------------

GEN_DATA_PATH="$(pwd)/../bins/gen_data"

if [ ! -x "$GEN_DATA_PATH" ]; then
    echo "Erro: executável $GEN_DATA_PATH não encontrado ou não tem permissão de execução."
    exit 1
fi

# -------------------------------
# VERIFICAÇÃO DE PERMISSÕES
# -------------------------------

check_and_fix_permissions() {
    local dir="$1"
    if [ ! -w "$dir" ]; then
        echo "Permissão insuficiente para escrever em $dir."
        echo "Tentando corrigir permissões..."
        chmod -R u+rwX,go+rwX "$dir" 2>/dev/null

        if [ ! -w "$dir" ]; then
            echo "Erro: ainda sem permissão de escrita em $dir após tentativa de correção."
            echo "Verifique as permissões manualmente e execute o script novamente."
            exit 1
        else
            echo "Permissões corrigidas com sucesso em $dir."
        fi
    fi
}

check_and_fix_permissions "$DATA_DIR"
check_and_fix_permissions "$LOG_DIR"

# -------------------------------
# GERAR SEEDS ÚNICAS
# -------------------------------

generate_unique_seeds() {
    local total=$1
    local max_seed=1000000
    shuf -i 1-$max_seed -n $total
}

TOTAL_RUNS=$(( ${#N_SPINS[@]} * N_RUNS ))
SEEDS=($(generate_unique_seeds $TOTAL_RUNS))

# -------------------------------
# CRIAR ARQUIVO DE ARGUMENTOS
# -------------------------------

index=0
for N in "${N_SPINS[@]}"; do
    for ((j=1; j<=N_RUNS; j++)); do
        SEED=${SEEDS[$index]}
        echo "$N 30000 metropolis $SEED on" >> "$ARGS_FILE"
        ((index++))
    done
done

echo "Arquivo de argumentos criado em $ARGS_FILE com $TOTAL_RUNS linhas."

# -------------------------------
# FUNÇÃO PARA RODAR CADA JOB
# -------------------------------

run_job() {
    N=$1
    M=$2
    METHOD=$3
    SEED=$4
    TEST=$5

    LOGFILE="$LOG_DIR/log_N${N}_seed${SEED}.txt"

    echo "Executando: $GEN_DATA_PATH $N $M $METHOD $SEED $TEST" > "$LOGFILE"

    ATTEMPTS=0
    MAX_ATTEMPTS=3

    while (( ATTEMPTS < MAX_ATTEMPTS )); do
        "$GEN_DATA_PATH" "$N" "$M" "$METHOD" "$SEED" "$TEST" >> "$LOGFILE" 2>&1
        EXIT_CODE=$?
        if [ $EXIT_CODE -eq 0 ]; then
            echo "Finalizado com sucesso." >> "$LOGFILE"
            break
        else
            echo "Erro detectado (tentativa $((ATTEMPTS+1)) ). Tentando novamente..." >> "$LOGFILE"
            ((ATTEMPTS++))
            sleep 2
        fi
    done

    if (( ATTEMPTS == MAX_ATTEMPTS )); then
        echo "Falhou após $MAX_ATTEMPTS tentativas." >> "$LOGFILE"
    fi
}

export -f run_job
export LOG_DIR
export GEN_DATA_PATH

# -------------------------------
# EXECUTAR EM PARALELO
# -------------------------------

NUM_CORES=$(nproc --all)

cat "$ARGS_FILE" | parallel -j "$NUM_CORES" --colsep ' ' run_job {1} {2} {3} {4} {5}

echo "Todas as tarefas enviadas."
