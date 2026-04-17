#!/usr/bin/env bash

# on a mac system this script can be sym-linked into /usr/local/bin
# giving the ability to launch the the shiny app by typing
# `chimera` at any command prompt regardless of directory.
# this needs to be in the top level of the shiny app directory
# to function that way.

#Get source dir

SOURCE="${BASH_SOURCE[0]}"
while [ -L "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"

echo $SCRIPT_DIR

R -e "shiny::runApp('${SCRIPT_DIR}/app.R', launch.browser=TRUE)"