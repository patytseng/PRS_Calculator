#####======================================== Imports ========================================#####

#built-in
import os
import gc
import sys
import gzip
import time
import math
import numpy as np
import argparse
import concurrent.futures

#installed
import allel
import pandas as pd
from psutil import virtual_memory

#####======================================== Input Arguments ========================================#####

#arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Lets get some PRS done~')
    parser.add_argument('--vcf', type=str, nargs='+', required=True, help='file path to directory containing vcf.gz files (can be multiple but number of files should be some multiple of 22)')
    parser.add_argument('--suffix', type=str, required=True, help='common suffix of all vcf files')
    parser.add_argument('--PRS', type=str, required=True, help='file path to polygenic risk scoring file (takes a csv with ONE header line)')
    parser.add_argument('--CHR', type=str, required=True, help='column header of chromosome number')
    parser.add_argument('--BP', type=str, required=True, help='column header of base pair')
    parser.add_argument('--EFFECT', type=str, required=True, help='column header of effect allele')
    parser.add_argument('--REFERENCE', type=str, required=True, help='column header of reference allele')
    parser.add_argument('--WEIGHT', type=str, required=True, help='column header of effect size or weight of each SNP')
    parser.add_argument('--ncore', type=int, default=os.cpu_count(), help='number of cores (will use all if not specified)')
    parser.add_argument('--nthread', type=int, default=1, help='number of threads per parallel process (increase if CPU utilization is low)')
    parser.add_argument('--out', type=str, required=False, help='output file name (will be a .csv file)')
    parser.add_argument('--explicit', default=False,type=bool, help='long logging and debug information')
    return parser.parse_args()

def display_args(args):
    variables = vars(args)
    print('\n\n\nInput Options:')
    for argument in variables:
        print(f'  --{argument} {variables[argument]}')

#####======================================== Grabbing VCF File Paths ========================================#####
    
# grabbing chromosome identifier from file path
def chromosome_label(path,prefix):
    path = list((path.split('/')[-1]).split('.')[:-2])
    index = len(prefix)
    for item in path:
        if item[0:index] == prefix:
            if item[len(prefix):].isnumeric():
                return (int(item[index:]))
            else:
                return None

# generate list of all compressed vcf files
def get_vcf_file_path_list(vcf_path,suffix):
    output = []
    for i in range(0,22):
        output.append([])
    
    vcf_file_path_list = []
    # grabbing all compressed .vcf.gz files and their directory
    for path in vcf_path:
        for filename in os.listdir(path):
            if filename.endswith(suffix) and path[-1] == "/":
                vcf_file_path_list.append(path + filename)
            elif filename.endswith(suffix):
                vcf_file_path_list.append(path + "/" + filename)
    
    # identifying chromosome number prefix
    vcf_file_path_list.sort()
    mask = (vcf_file_path_list[0].split('/')[-1]).split('.')[:-2]
    for path in vcf_file_path_list[1:22]:
        path = (path.split('/')[-1]).split('.')[:-2]
        for index in range(0,len(mask)):
            if mask[index] == path[index]:
                mask[index] = 0
    
    # removing common elements
    removal_list = []
    for item in mask:
        if item == 0:
            removal_list.append(item)
    for item in removal_list:
        mask.remove(item)
    
    prefix = mask[0][:-1]
    
    # removing non-chr#s 1-22
    removal_list = []
    
    for path in vcf_file_path_list:
        if chromosome_label(path,prefix) == None:
            removal_list.append(path)
    for path in removal_list:
        vcf_file_path_list.remove(path)
        
    # sorting all obtained paths
    vcf_file_path_list.sort()
    
    # grouping obtained vcf file paths
    for path in vcf_file_path_list:
        output[chromosome_label(path,prefix)-1].append(path)
    return output

#####======================================== PRS Grabber ========================================#####

def get_PRS_targets(path,CHR,BP,EFFECT,REFERENCE,WEIGHT):

    # generating a blank list
    PRS_targets = []

    # opening file
    df = pd.read_csv(path)
    
    #reading into array
    for index, row in df.iterrows():
        PRS_targets.append((row[CHR], row[BP],row[EFFECT],row[REFERENCE],row[WEIGHT]))
    
    return PRS_targets

#####======================================== Subject ID List ========================================#####

# obtaining subject ID list
def get_ID_library(vcf_file_path):
    ID_library = {}
    
    for path in vcf_file_path:
        # pulling relevant line from vcf file header
        with gzip.open(path,'rt') as file:
            for line in file:
                line = list(line[:-1].split('\t'))
                if line[0] == '#CHROM':
                    index = 1
                    for item in line:
                        if item == 'FORMAT':
                            break
                        else:
                            index += 1
                    IDs = line[index:]
                    break

        # converting to dictionary for easier processing
        for ID in IDs:
            ID_library[ID] = float(0)

    return ID_library

#####======================================== Parallel Tasks ========================================#####

class get_vcf_line(object):
    def __init__(self, vcf_file_path, chromosome_number, position):
        self.full_data = allel.read_vcf(vcf_file_path,fields=['variants/ALT','variants/REF','calldata/GT','samples'],region=f'{chromosome_number}:{position}-{position}')
    def __enter__(self):
        return self.full_data
    def __exit__(self ,type, value, traceback):
        if type == dict:
            self.full_data.clear()
        del self.full_data
        gc.collect()

def chunk_PRS_targets(PRS_targets,chunksize):
    
    # generating empty lists
    PRS_target_chunks = []
    chunk = []
    
    # chunking the lists
    loopcount = 1
    for target in PRS_targets:
        chunk.append(target)
        if loopcount % chunksize == 0:
            PRS_target_chunks.append(chunk)
            chunk = []
        loopcount += 1
        
    # cleaning up last chunk
    PRS_target_chunks.append(chunk)
    
    return PRS_target_chunks

def get_partial_PRS(vcf_file_paths,PRS_targets,ID_library,explicit):
    # tracking SNP casess
    cases = [[0 for foo in range(0,22)] for foo in range(0,4)]
    
    # iterating over portion of PRS
    for SNP_info in PRS_targets:
        # grabbing line from vcf file
        for path in vcf_file_paths[int(SNP_info[0])-1]:
            with get_vcf_line(path,SNP_info[0],SNP_info[1]) as line:
                if line == None:
                    cases[0][SNP_info[0]-1] += 1
                    if explicit:
                        print(f'position {SNP_info[1]} on CHR {SNP_info[0]} is missing in\n{path}')
                elif line['variants/ALT'][0][0] == SNP_info[2] and line['variants/REF'][0] == SNP_info[3]:
                    cases[1][SNP_info[0]-1] += 1
                    
                    # summing in case of correct effect/reference
                    index = 0
                    for ID in line['samples']:
                        ID_library[ID] = float(ID_library[ID] + float(SNP_info[4]) * sum(line['calldata/GT'][0][index]))
                        index += 1
                    
                elif line['variants/ALT'][0][0] == SNP_info[3] and line['variants/REF'][0] == SNP_info[2]:
                    cases[2][SNP_info[0]-1] += 1
                    if explicit:
                        print(f'position {SNP_info[1]} on CHR {SNP_info[0]} has reversed alleles in\n{path}')
                    
                    # summing in case of reversed effect/reference
                    index = 0
                    for ID in line['samples']:
                        ID_library[ID] = float(ID_library[ID] + float(SNP_info[4]) * (2-sum(line['calldata/GT'][0][index])))
                        index += 1
                    
                else:
                    cases[3][SNP_info[0]-1] += 1
                    if explicit:
                        print(f'position {SNP_info[1]} on CHR {SNP_info[0]} has mismatched alleles in\n{path}')
                    
    # accounting for multi-file systems
    file_num = len(vcf_file_paths[0])
    if ( sum( sum(case%file_num for case in cases[index]) for index in range(0,4)) == 0):
        for index in range(0,4):
            cases[index] = [int(case/file_num) for case in cases[index]]
    
    return [ID_library,cases]
    
def task(vcf_file_path,PRS_targets,ID_library,chunk_number,total_chunks,threads,explicit):
    
    #print(f'Chunk {chunk_number} of {total_chunks} started')
    
    cases = [[0 for foo in range(0,22)] for foo in range(0,4)]

    PRS_target_chunks = chunk_PRS_targets(PRS_targets,threads)
                    
    for PRS_targets in PRS_target_chunks:
    
        # launching threads
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            future_list = []
            chunk_counter = 1
            total_chunks = len(PRS_targets)
            
            for PRS_target in PRS_targets:
                future_list.append(executor.submit(get_partial_PRS,vcf_file_path,[PRS_target],ID_library,explicit))
            
            chunk_counter += 1
    
        # extracting threaded calculation results and summing them
        for future in future_list:
            for key in ID_library:
                ID_library[key] = ID_library[key] + future.result()[0][key]
        
        # extracting threaded calculation cases
        for future in future_list:
            for x in range(0,4):
                for y in range(0,22):
                    cases[x][y] += future.result()[1][x][y]
        
        # cleanup
        future_list = None
        gc.collect()
        
    #print(f'Chunk {chunk_number} of {total_chunks} completed')
    
    return [ID_library,cases]
    

#####======================================== Output ========================================#####

def get_output_name(prefix,suffix,vcf_file_paths,PRS_file_path):
    # convert paths to lists of strings of file name
    VCF_identifier = list(vcf_file_paths[0][len(prefix):-1*len(suffix)].split('.'))
    # identifies incongruities between file names and replaces with None
    if len(vcf_file_paths) > 1:
        for path in vcf_file_paths[1:]:
            split_path = list(path[len(prefix):-1*len(suffix)].split('.'))
            for index in range(0,len(VCF_identifier)):
                # removes items that are not in common between names
                if VCF_identifier[index] != split_path[index]:
                    VCF_identifier[index] = 0
    
    # removing None items and chr#
    remove_list = []
    for item in VCF_identifier:
        if item == 0:
            remove_list.append(item)
        elif item[0:3] == 'chr':
            remove_list.append(item)
    for item in remove_list:
        VCF_identifier.remove(item)
    
    VCF_identifier = '.'.join(VCF_identifier)
    
    # final checks just in case...
    if VCF_identifier[0] == '/':
        VCF_identifier = VCF_identifier[1:]
    
    PRS_identifier = '.' + '.'.join(list((PRS_file_path.split('/')[-1]).split('.')[:-1])) + '_score'
        
    return VCF_identifier + PRS_identifier + '.csv'

def write_output(ID_library,out_name):
    
    # opening file to write to
    try:
        output_file = open(out_name,'x')
    except:
        output_file = open(out_name,'w')
    
    # writing to file: header then data
    output_file.write('ID,Score\n')
    for key in ID_library:
        output_file.write(f'{key},{ID_library[key]}\n')
    
#------------------------------------------------------------------------------------------------------#
#####======================================== Main Program ========================================#####
#------------------------------------------------------------------------------------------------------#

def main():
    
#======================================== Taking Inputs (Preparation)
    
    # setting up input arguments
    args = parse_args()
    display_args(args)
    
    # initiating timer
    tic = time.time()
    
    # grabbing vcf file names in vcf file directory
    vcf_file_path_list = get_vcf_file_path_list(args.vcf,args.suffix)
    
    # displaying list of vcf files
    print('\n\n\nObtained list of vcf files:')
    for paths in vcf_file_path_list:
        for path in paths:
            print(path)
    
    # processing PRS file to obtain list of SNPs to target
    raw_PRS_targets = get_PRS_targets(args.PRS,args.CHR,args.BP,args.EFFECT,args.REFERENCE,args.WEIGHT)
    SNP_count = len(raw_PRS_targets)
    SNP_size = sys.getsizeof(raw_PRS_targets)
    print(f'\n\n\nGathered [{SNP_count} total SNPs from PRS file]')
    print(f'SNP list uses {SNP_size} Bytes')
    
    # obtaining list of subject IDs
    ID_library = get_ID_library(vcf_file_path_list[0])
    ID_count = len(ID_library)
    ID_size = sys.getsizeof(ID_library)
    print(f'\nGathered {ID_count} Subject IDs into ID library')
    print(f'ID library uses {ID_size} Bytes')
            
    # generating chunks of PRS_targets
    PRS_target_chunks = chunk_PRS_targets(raw_PRS_targets,math.ceil(len(raw_PRS_targets)/args.ncore))
    print('\nPRS targets Chunked for parallel processing')
    
    del raw_PRS_targets

    # time check
    toc = time.time()
    print(f"\nPreparations completed in {toc - tic:0.4f} seconds")

#======================================== PRS Calculation
    
    # multithreading over chunks of PRS targets
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncore) as executor:
        future_list = []
        chunk_counter = 1
        total_chunks = len(PRS_target_chunks)
        
        for PRS_targets in PRS_target_chunks:
            future_list.append( executor.submit( task, vcf_file_path_list, PRS_targets, ID_library, chunk_counter, total_chunks,  args.nthread, args.explicit ) )
            
            chunk_counter += 1
    
    toc = time.time()
    print(f'\nValue retrieval from VCF completed in {toc - tic:0.4f} seconds')

#======================================== Merging Multiple Outputs

    # extracting parallel calculation results and summing them
    PRS_ID_library_list = []
    for future in future_list:
        for key in ID_library:
            ID_library[key] = ID_library[key] + future.result()[0][key]
            
    # extracting parallel calculation cases
    cases = [[0 for foo in range(0,22)] for foo in range(0,4)]
    
    for future in future_list:
        for x in range(0,4):
            for y in range(0,22):
                cases[x][y] += future.result()[1][x][y]
    
#======================================== Displaying & Generating Results
    
    # displaying number of cases separated by chromosome:
    print('\nThese issues came up with the SNPs or Alleles:')
    for chromosome in range(1,23):
        print(f'\nChr{chromosome}:{" "*(3-len(str(chromosome)))}{cases[0][chromosome-1]} missing SNPs\n{" "*7}{cases[1][chromosome-1]} matched alleles\n{" "*7}{cases[2][chromosome-1]} reversed alleles\n{" "*7}{cases[3][chromosome-1]} mismatched alleles')
    
    # totals
    print(f'\nTotal: {sum(cases[0])} missing SNPs\n{" "*7}{sum(cases[1])} matched alleles\n{" "*7}{sum(cases[2])} reversed alleles\n{" "*7}{sum(cases[3])} mismatched alleles')
    
    # summary of used & excluded
    snps_used = (sum(cases[1])+sum(cases[2]))
    snps_excluded = (sum(cases[0])+sum(cases[3]))
    print(f'\nUsed {round(100*snps_used/SNP_count,3)}% of the SNPs ({snps_used}/{SNP_count})')
    print(f'\nExcluded {round(100*snps_excluded/SNP_count,3)}% of the SNPs ({snps_excluded}/{SNP_count})')
    
    # generating output file name
    if args.out == None:
        out_name = get_output_name(args.vcf[0],'.vcf.gz',vcf_file_path_list[0],args.PRS)
    else:
        out_name = args.out + '.score.csv'
    
    # writing output
    print('\nWriting output...')
    write_output(ID_library,out_name)
    print(f'Results written to {out_name}\n')

    # closing timer and collecting result
    toc = time.time()
    print(f"Completed PRS Calculation in {toc - tic:0.4f} seconds")

#####======================================== Initiation ========================================#####
    
if __name__ == '__main__':
    # initiating main program
    main()