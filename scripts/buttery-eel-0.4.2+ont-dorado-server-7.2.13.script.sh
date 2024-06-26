#!/bin/bash

PORT=4000

# https://unix.stackexchange.com/questions/55913/whats-the-easiest-way-to-find-an-unused-local-port
get_free_port1(){
    PORT=$(netstat -aln | awk '
    $6 == "LISTEN" {
        if ($4 ~ "[.:][0-9]+$") {
        split($4, a, /[:.]/);
        port = a[length(a)];
        p[port] = 1
        }
    }
    END {
        for (i = 5000; i < 65000 && p[i]; i++){};
        if (i == 65000) {exit 1};
        print i
    }')
}

get_free_port2() {
	for port in $(seq 5000 65000); do
		echo "trying port $port" >&2
		PORT=$port
		ss -lpna | grep -q ":$port " || break
	done
}

get_free_port(){
    if which netstat > /dev/null 2> /dev/null;
    then
        get_free_port1
    else
        get_free_port2
    fi
    test -z "${PORT}" && PORT=60000
    echo "Using port ${PORT}"
}

REALPATH=$(dirname "$(readlink -f "$0")")
test -z ${REALPATH} && echo "REALPATH could not be deduced" && exit 1
source ${REALPATH}/../venv3/bin/activate || { echo "${REALPATH}/../venv3 could not be activated" && exit 1; }
get_free_port
buttery-eel "$@" -g /data/install/buttery-eel-0.4.2+dorado7.2.13/ont-dorado-server-7.2.13/bin --port ${PORT} --use_tcp
exit_code=$?
deactivate
exit "$exit_code"
