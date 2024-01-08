import os
import subprocess
import time
import json
import csv
import time
import sys
import multiprocessing
from multiprocessing import Pool
from os import path

DIR = os.getcwd() + '/Desktop/chloralign/'
SRC = DIR + 'source/'
TMP = DIR + 'temp/'
CSV = DIR + 'AccessionList2.txt'
FASTP = SRC + 'fastp'
REFGENOME = 'HA89_cytoplasm.fasta'
OUT = DIR + 'output/'
ERR = 'errorReport.txt'
numProcessors = 50
maxTasksPerChild = 1
argList = []

def command(cmd: str):
 subprocess.run([cmd], shell = True)

def timestamp():
 return str('[' + time.ctime() + ']: ')

def argumentCheck(argList):
 for arg in argList:
  if argList[0] == '-dir' or argList[0] == '-d':
   DIR = argList[1]
   SRC = DIR + 'source/'
   OUT = DIR + 'output/'
   TMP = DIR + 'temp/'
   CSV = SRC + 'AccessionList2.txt'
   FASTP = SRC + 'fastp'
  elif argList[0] == '-referencegenome' or argList[0] == '-r':
   REFGENOME = argList[1]
  elif argList[0] == '-input' or argList[0] == '-i':
   CSV = SRC + argList[1]
  elif argList[0] == '-numprocessors' or argList[0] == '-np':
   numProcessors = int(argList[1])
  del argList[0]
  del argList[0]

def processRun(run):
 try:
  with open (OUT + run + '.lite.sam') as currentRunFile:
   print(timestamp() + run + '.lite.sam found.')
 except IOError:
  print(timestamp() + run + ' download starting.')
  filename = str(run) + 'lite.1'
  filepath = TMP + str(filename)
  command(SRC + 'sratoolkit.3.0.2-centos_linux64/bin/prefetch ' + run + ' -O ' + TMP + run + ' --max-size 420000000000 2> ' + ERR + ' > dump.txt')
  command(SRC + 'sratoolkit.3.0.2-centos_linux64/bin/parallel-fastq-dump-0.6.5 --split-files ' + run + ' -O ' + TMP + ' -t ' + TMP + ' -f 2> ' + ERR )
  command(FASTP + ' -h ' + TMP + ' -j ' + TMP + ' -i ' + TMP + run + '_1.fastq -I ' + TMP + run + '_2.fastq -o ' + TMP + run + '_1.trim.fastq -O ' + TMP + run + '_2.trim.fastq 2> ' + ERR)
  command('bwa mem -t 1 ' + SRC + REFGENOME + ' ' + TMP + run + '_1.fastq -I ' + TMP + run + '_2.trim.fastq > ' + OUT + run + '.lite.sam 2> ' + ERR)
  command('samtools view -b -F 1 ' + TMP + run + '.lite.sam -o ' + OUT + run + '.lite.sam' + ' 2> ' + ERR)
  command('rm -rf ' + TMP + run + '*.fastq ' + TMP + run)
 print(timestamp() + run + ' download complete.')

def processBiosample(biosample):
  command('samtools sort -T ' + TMP + ' ' + TMP + biosample + '.bam' + ' -o ' + TMP + biosample + '.bam')
  print(biosample + '.bam file created.')
  command(SRC + 'bcftools-1.18/bcftools mpileup -O u -f ' + SRC + REFGENOME + ' ' + TMP + biosample + '.bam | ' + SRC + 'bcftools-1.18/bcftools call -c -O v --ploidy 1 | bcftools view -i \'INFO/DP>10\' -Oz -o ' + OUT + biosample + '.vcf.gz')
  command(SRC + 'bcftools-1.18/bcftools index -t -f ' + TMP + biosample + '.vcf.gz')
  print(timestamp() + biosample + ' load complete.')


def list_homology(run_list, biosample_list, run_sample_dict):
 print(timestamp() + 'Merging runs by biosample ID...')
 biosample_dict = {}
 homologous_sample_list = []
 for biosample in biosample_list:
  for run in run_sample_dict:
   if (run_sample_dict[run] == biosample):
    homologous_sample_list.append(run)
  biosample_dict[biosample] = homologous_sample_list
  homologous_sample_list = []
 return biosample_dict

def check_AccessionList():
 try:
  with open(CSV) as f:
   print(timestamp() + CSV + ' located.')
 except IOError:
  print(timestamp() + 'ERROR: Input file not detected. Check source directory.')
  command('ls')
  print('couldnt find ' + CSV)
  quit()

def check_Fastp():
 try:
  with open(FASTP) as f:
   print(timestamp() + 'Fastp file located.')
 except IOError:
  print(timestamp() + 'Fastp file not detected. Downloading to source directory...')
  command('wget -P ' + SRC + ' http://opengene.org/fastp/fastp')
  print(timestamp() + 'Success: Fastp file downloaded to source directory.')

def check_SRAToolkit():
 try:
  with open(SRC + 'sratoolkit.3.0.2-centos_linux64.tar.gz') as f:
   print(timestamp() + 'SraToolkit file located.')
 except IOError:
  print(timestamp() + 'SraToolkit not detected. Downloading to source directory...')
  command('wget -P ' + SRC + ' https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-centos_linux64.tar.gz')
  print(timestamp() + 'Success: Sratoolkit downloaded.')
  print(timestamp() + 'Decompressing SraToolkit...')
  command('tar -xzf sratoolkit.3.0.2-centos_linux64.tar.gz')
  print(timestamp() + 'Success: Sratoolkit ready.')

def main():
 print('\n=====\nWelcome to Chloralign v2.0\n=====\n')
 time_start = time.time()
 for arg in sys.argv:
  argList.append(arg)
 del argList[0]
 argumentCheck(argList)
 print(timestamp() + 'Checking for mandatory modules...')
 check_AccessionList()
 check_Fastp()
 check_SRAToolkit()
 with open(CSV, newline = '', encoding='utf-8-sig') as filename:
  print(timestamp() + 'Input file opened.')
  reader = csv.DictReader(filename)
  file = csv.reader(filename)
  run_list = []
  biosample_list = []
  next(filename)
  print(timestamp() + 'Filtering out samples without paired reads.')
  print(timestamp() + 'Filtering out samples collected using PACBIO.')
  for row in file:
   row = str(row).split(',')
   if row[11] != 'PACBIO_SMRT' and row[12] != 'PACBIO_SMRT' and row[7] != 'SINGLE':
    run = str(row[0].lstrip('[\''))
    run_list.append(run)
    biosample = row[2]
    biosample_list.append(biosample)
 run_sample_dict = dict(zip(run_list, biosample_list))
 biosample_dict = list_homology(run_list, biosample_list, run_sample_dict)
 command('chmod a+x ' + SRC + 'fastp')
 command('bwa index ' + SRC + REFGENOME + ' 2> ' + ERR)
 print(timestamp() + 'Initialization complete.')
 print(timestamp() + 'Starting run processing...')
 pool = multiprocessing.Pool(processes=numProcessors)
 pool = multiprocessing.pool.Pool(maxtasksperchild=maxTasksPerChild)
 print(timestamp() + 'There are currently 1 processors in use.')
 with Pool(numProcessors, maxtasksperchild = maxTasksPerChild) as pool:
  pool.map(processRun, run_list)
 print(timestamp() + 'Confirmed presence of necessary .lite.sam files.')
 command('rm -rf ' + TMP + '*')
 for biosample in biosample_list:
  homologous_run_list = biosample_dict[biosample]
  if len(homologous_run_list) > 1:
   for homo_run in homologous_run_list:
    command('samtools view -b -F 1 ' + OUT + run + '.lite.sam -o ' + TMP + run + '.bam')
   command('samtools merge ' + biosample_bam + (' ' + TMP).join(homologous_run_list) + ' -f')
   print(timestamp() + 'Identical biosample ID ' + biosample + ' detected: Data merged.')
  else:
   command('samtools view -b -F 1 ' + OUT + run + '.lite.sam -o ' + TMP + biosample + '.bam')
 with Pool(numProcessors, maxtasksperchild = maxTasksPerChild) as pool:
  pool.map(processBiosample, biosample_list)
 command(SRC + 'bcftools-1.18/bcftools merge -o ' + OUT + 'data.vcf.gz ' + TMP + ('.vcf.gz ' + TMP).join(biosample_list) + '.vcf.gz')
 time_end = time.time()
 time_total = time_end - time_start
 print(timestamp() + 'Program complete. (Time = ' + str(time_total) + ' minutes)')

main()
