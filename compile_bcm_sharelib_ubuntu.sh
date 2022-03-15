#!/bin/bash
echo 'Building BCM library...'


g++ -Wall -fPIC -shared bcm.cpp -o libbcm.o -std=c++0x



