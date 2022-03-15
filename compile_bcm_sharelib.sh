#!/bin/bash

echo 'Building BCM library...'
gcc -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0,-install_name,@loader_path/libbcm.dylib -o libbcm.dylib bcm.cpp



