#!/bin/bash
sleep $((3 * 3600))  # 3 hours
./scripts/submit_scripts/submit_pr.sh -s luc -e fpichard@umn.edu -S 00 --n-subjects 182 -c pr_genparams,pr_simulate,pr_recovery -k igt -m vppdecay
./scripts/submit_scripts/submit_pr.sh -s luc -e fpichard@umn.edu -S 00 --n-subjects 182 -c pr_genparams,pr_simulate,pr_recovery -k igt -m vppboth
