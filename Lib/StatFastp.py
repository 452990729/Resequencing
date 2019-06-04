#!/usr/bin/env python2


import sys
import re
import os
import json
from glob import glob

def ExtractJson(file_in):
    with open(file_in, 'r') as f:
        dict_json = json.load(f)
    before =  dict_json['summary']['before_filtering']
    after = dict_json['summary']['after_filtering']
    before_total = before['total_reads']
    after_total = after['total_reads']
    pecent = round(float(after_total)*100/before_total, 2)
    berfore_q20 = before['q20_rate']
    berfore_q30 = before['q30_rate']
    after_q20 = after['q20_rate']
    after_q30 = after['q30_rate']
    return '\t'.join([str(before_total), str(after_total), str(pecent),\
                      str(berfore_q20), str(berfore_q30),\
                      str(after_q20), str(after_q30)])

def main():
    print 'before_total\tafter_total\tfilter(%)\tberfore_q20\tberfore_q30\tafter_q20\tafter_q30'
    for fl in sorted(glob(sys.argv[1]+'/*_QC_report.json')):
        lb = os.path.basename(fl).rstrip('_QC_report.json')
        print lb+'\t'+ExtractJson(fl)


if __name__ == '__main__':
    main()
