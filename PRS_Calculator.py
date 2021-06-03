#####======================================== Imports ========================================#####

import os
import gzip
import allel
import concurrent.futures
import time
import math
import argparse
import pandas as pd

#####======================================== Input Check ========================================#####

def display_args(args):
    print(f'\n\n\nInput Options:\n--vcf {args.vcf}\n--suffix {args.suffix}\n--PRS {args.PRS}\n--CHR {args.CHR}\n--BP {args.BP}\n--EFFECT {args.EFFECT}\n--REFERENCE {args.REFERENCE}\n--WEIGHT {args.WEIGHT}\n--ncore {args.ncore}\n--out {args.out}\n\n\n')

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
            ID_library[ID] = 0

    return ID_library

#####======================================== Parallel Tasks ========================================#####

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

def get_vcf_line(vcf_file_path,chromosome_number,position):
    return allel.read_vcf(vcf_file_path,fields=['variants/ALT','variants/REF','calldata/GT','samples'],region=f'{chromosome_number}:{position}-{position}')

def get_partial_PRS(vcf_file_paths,PRS_targets,ID_library):
        
    # iterating over portion of PRS
    for SNP_info in PRS_targets:
        
        # grabbing line from vcf file
        for path in vcf_file_paths[int(SNP_info[0])-1]:
            line = get_vcf_line(path,SNP_info[0],SNP_info[1])
            if line == None:
                print(f'position {SNP_info[1]} on CHR {SNP_info[0]} not found in\n{path}')
            elif line['variants/ALT'][0][0] == SNP_info[2] and line['variants/REF'][0] == SNP_info[3]:
                index = 0
                for ID in line['samples']:
                    ID_library[ID] = ID_library[ID] + float(SNP_info[4]) * sum(line['calldata/GT'][0][index])
                    index += 1
            else:
                print(f'position {SNP_info[1]} on CHR {SNP_info[0]} does not have matching alleles in\n{path}')
                
    return ID_library
    
def task(vcf_file_path,PRS_targets,ID_library,chunk_number,total_chunks):
    print(f'Chunk {chunk_number} of {total_chunks} started')
    output = get_partial_PRS(vcf_file_path,PRS_targets,ID_library)
    print(f'Chunk {chunk_number} of {total_chunks} completed\n')
    return output

#####======================================== Output ========================================#####

def get_output_name(prefix,suffix,vcf_file_paths):
    # convert paths to lists of strings of file name
    output = list(vcf_file_paths[0][len(prefix):-1*len(suffix)].split('.'))
    # identifies incongruities between file names and replaces with None
    if len(vcf_file_paths) > 1:
        for path in vcf_file_paths:
            split_path = list(path[len(prefix):-1*len(suffix)].split('.'))
            for index in range(0,len(output)):
                # removes items that are not in common between names
                if output[index] != split_path[index]:
                    output[index] = 0
    
    # removing None items and chr#
    remove_list = []
    for item in output:
        if item == 0:
            remove_list.append(item)
        elif item[0:3] == 'chr':
            remove_list.append(item)
    for item in remove_list:
        output.remove(item)
    
    output = '.'.join(output)
    
    # final checks just in case...
    if output[0] == '/':
        output = output[1:]
    return output + '.scores.csv'

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

#####======================================== Main Program ========================================#####
        
def main():
    # setting up argument parsing
    parser = argparse.ArgumentParser(description='Lets get some PRS done~ Note: VCFs are expected to be gzipped')
    parser.add_argument('--vcf', type=str, nargs='+', help='file path to vcf file directory')
    parser.add_argument('--suffix', type=str, help='common suffix of all vcf files')
    parser.add_argument('--PRS', type=str, help='file path to polygenic risk scoring file')
    parser.add_argument('--CHR', type=str, help='column header of chromosome number')
    parser.add_argument('--BP', type=str, help='column header of base pair')
    parser.add_argument('--EFFECT', type=str, help='column header of effect allele')
    parser.add_argument('--REFERENCE', type=str, help='column header of reference allele')
    parser.add_argument('--WEIGHT', type=str, help='column header of effect size or weight of each SNP')
    parser.add_argument('--ncore', type=int, help='number of cores')
    parser.add_argument('--out', type=str, help='output file name (will be a .csv file)')
    args = parser.parse_args()
    
    display_args(args)
    
    # initiating timer
    tic = time.time()
    
    # collecting basic variables
    vcf_file_path_list = get_vcf_file_path_list(args.vcf,args.suffix)
    
    print('Obtained list of vcf files:')
    for paths in vcf_file_path_list:
        for path in paths:
            print(path)
    print('\n\n')
    
    raw_PRS_targets = get_PRS_targets(args.PRS,args.CHR,args.BP,args.EFFECT,args.REFERENCE,args.WEIGHT)
    print('PRS target SNP information gathered')
    
    ID_library = get_ID_library(vcf_file_path_list[0])
    print('Subject IDs Gathered')
    
    # naming output file
    if args.out == None:
        out_name = get_output_name(args.vcf[0],'.vcf.gz',vcf_file_path_list[0])
    else:
        out_name = args.out + '.score.csv'
            
    # generating chunks of PRS_targets
    PRS_target_chunks = chunk_PRS_targets(raw_PRS_targets,math.ceil(len(raw_PRS_targets)/args.ncore))
    print('PRS targets Chunked for parallel processing')
    
    # time check
    toc = time.time()
    print(f"\nPreparations completed in {toc - tic:0.4f} seconds\n")
    
    # multithreading over sections of PRS targets
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncore) as executor:
        future_list = []
        chunk_counter = 1
        total_chunks = len(PRS_target_chunks)
        
        for PRS_targets in PRS_target_chunks:
            future = executor.submit(task,vcf_file_path_list,PRS_targets,ID_library,chunk_counter,total_chunks)
            future_list.append(future)
            
            chunk_counter += 1
            
        PRS_ID_library_list = []
        for future in future_list:
            PRS_ID_library_list.append(future.result())
    
    # summing results of each calculated PRS chunk
    for PRS_ID_library in PRS_ID_library_list:
        for key in ID_library:
            ID_library[key] = round(ID_library[key] + PRS_ID_library[key],4)
        
    # writing output
    print('\nWriting output...')
    write_output(ID_library,out_name)
    print(f'Results written to {out_name}\n')

    # closing timer and collecting result
    toc = time.time()
    print(f"Completed PRS Calculation in {toc - tic:0.4f} seconds")

if __name__ == '__main__':
    # initiating main program
    main()