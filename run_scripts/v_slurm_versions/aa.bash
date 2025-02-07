#!/bin/bash

callingDir="$(pwd)"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

sd=$(dirname "$0")
sd2=$(dirname $(readlink -f $0))

echo "$callingDir"
echo "$SCRIPT_DIR"
echo "$sd"
echo "$sd2"
