#!/usr/bin/python

import random
lines = open('greencap.dat').readlines()
random.shuffle(lines)
open('greencap.dat', 'w').writelines(lines)
