#!/usr/bin/env python3
import sys
import statistics

lines = sys.stdin.read().splitlines()
chroms = [x.split()[0] for x in lines]

distance_estimations = [float(x.split()[1]) for x in lines[:-1]]
positions = [int(x.strip('chrA1.')) for x in chroms]
distances = [a - b for a, b in zip(positions, positions[1:])]

assert len(distances) == len(distance_estimations)
ks = [int(d / est) for (d, est) in zip(distances, distance_estimations)
      if est > 0]


print(statistics.pvariance(ks))
