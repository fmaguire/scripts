#!/bin/sh

SOURCE=$(dd if=/dev/urandom count=1 2>/dev/null | openssl dgst -md5)
MACADDR="$(echo ${SOURCE} | \
sed 's/^\(..\)\(..\)\(..\)\(..\)\(..\)\(..\).*$/\1:\2:\3:\4:\5:\6/')"

for i in $(ifconfig -s | sed 1d - | cut -d ' ' -f 1); do
    echo "Downing interface '$i'... "
    ifconfig $i down
    echo "Changing MAC of '$i' to ${MACADDR}... "
    macchanger -m "$MACADDR" "$i"
    echo "Raising interface '$i'... "
    ifconfig $i up
done

