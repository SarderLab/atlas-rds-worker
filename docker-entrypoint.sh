#!/usr/bin/env bash

if [ "$1" = "--list_cli" ]; then
    cat list_cli.json
else
    if [ "$2" = "--xml" ]; then
        cat "$1".xml
    else
        Rscript "$1".R "${@:2}"
    fi
fi
