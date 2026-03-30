#!/bin/bash

# This small script is intended to streamline developing analysis. The assumption here is
# that the code is developed in one location and then tested in another. Repeated testing
# requires copying the updated files to the testing location. The script simply provides
# the (otherwise trivial) copying formula so one never needs to type it in anew.
#
# NOTE: The address of the development repo needs to be hardcoded for any indivisual user.

DEVEL_REPO=~/SMASH-ROOT-analysis

# copy
cp "$DEVEL_REPO/src/"* ../src/ && cp "$DEVEL_REPO/input/"* ../input/

