---
name: Bug report
about: Create a report to help us improve
title: 'SnakeFileTest'
labels: ''
assignees: ''

---

**Describe the bug**
Unable to run last unmapped_pie rule of the snakemake file due to errors with packages imported in python script 'blastresults_piechart.py'

**To Reproduce**
Steps to reproduce the behavior:
1. Go to 'eclip-qc'
2. Click on 'SnakeFileTest'
3. Scroll down to 'rule unmapped_pie'
4. See error

**Expected behavior**
If the rule was to run as intended, a piechart for each sample should be produced in user's pieChart directory. 

**Screenshots**
<img width="998" alt="Screen Shot 2022-12-16 at 2 50 33 PM" src="https://user-images.githubusercontent.com/98048034/208203288-633707b0-c083-4918-9732-9b6318ec59ed.png">
<img width="988" alt="Screen Shot 2022-12-16 at 3 09 01 PM" src="https://user-images.githubusercontent.com/98048034/208203289-2d72c468-db81-4c0c-ad39-033c3fda5bce.png">
<img width="994" alt="Screen Shot 2022-12-16 at 3 09 15 PM" src="https://user-images.githubusercontent.com/98048034/208203290-4d92771c-57b5-4d69-9e3d-fc53d44f7bb4.png">

**Additional context**
Seems to be an issue with a lot of the imported packages such as matplotlib and xlmodict, matplotlib is already listed as a dependency in the python3.yaml file we have written.
