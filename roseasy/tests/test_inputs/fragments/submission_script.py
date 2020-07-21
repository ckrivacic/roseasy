#!/usr/bin/python
#$ -S /usr/bin/python
#$ -N fragment_generation
#$ -o /Users/codykrivacic/cas/roseasy/roseasy/tests/test_inputs/fragments
#$ -e /Users/codykrivacic/cas/roseasy/roseasy/tests/test_inputs/fragments
#$ -cwd
#$ -r y
#$ -t 1-1
#$ -l arch=linux-x64
#$ -l mem_free=40G
#$ -l scratch=10G
#$ -l h_rt=6:00:00
#$ -q long.q

Non