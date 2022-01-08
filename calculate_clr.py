#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 10:55:36 2021

@author: ceciliethystrup
"""

import pandas as pd
import matplotlib.pyplot as plt
import os, sys
import numpy as np
import gzip
import seaborn as sns
import argparse
import numpy as np
import math


################### Create parser #########################
path = '/Users/ceciliethystrup/OneDrive/Dokumenter/DTU/Speciale E21/Metagenomics data/second_run/'

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--mapstat", nargs="*", help=".mapstat file", required=True)
parser.add_argument("-r", "--res", nargs="*", help=".res file", required=True)
parser.add_argument("-o", "--output", help="Output directory")
parser.add_argument("-d", "--db", help="Database to use for length templates", required=True)
parser.add_argument("-t", "--tax", help="Determines at what taxonomy-level FPKM-values should be calculated. Accepted values are 'species' or 'genus'.", required=True)


args = parser.parse_args()


##################### Define functions #######################
    
def split_id(row):
    ''' Split ID to get a shorter ID in .mapstat and .res file '''
    long_id = row['# id']
    long_id = str(long_id)
    short_id = long_id.split(' ')[0]
    return short_id
    
def split_ids(df):
    ''' Applies split_id-function on pandas dataframe'''
    df['# shortID'] = df.apply(split_id, axis=1)

def sum_fragment_count(row, tax_list, tax_dict, flag):
    ''' Sums up fragment count on either genus or species level'''
    if flag == 'species':
        species = row['species']
        fragmentCount = row['fragmentCountAln']

        if species in tax_list:
            if not species == 'nan':
                if species not in tax_dict:
                    tax_dict[species] = fragmentCount
                else:
                    if not fragmentCount == 'nan':
                        summed_fragment_count = int(tax_dict[species]) + fragmentCount
                        tax_dict[species] = summed_fragment_count
    else:
        genus = row['genus']
        fragmentCount = row['fragmentCountAln']

        if genus in tax_list:
            if not genus == 'nan':
                if genus not in tax_dict:
                    tax_dict[genus] = fragmentCount
                else:
                    if not fragmentCount == 'nan':
                        summed_fragment_count = int(tax_dict[genus]) + fragmentCount
                        tax_dict[genus] = summed_fragment_count

def geometric_mean(iterable):
    a = np.log(iterable)
    val = np.exp(a.mean())
    return val

def calc_clr(value, geom_mean):
    val = int(value)/geom_mean
    clr_val = math.log(val, 10)
    return clr_val


###############################################################################################


if len(args.mapstat) != 0 and len(args.res) != 0:
    for mapstat_file, res_file in zip(args.mapstat, args.res):
        #reset variables
        species_dict = {}
        genus_dict = {}
        clr_dict = {}
        
        ####
        sys.stdout.write("######### Working on {} #########\n".format(mapstat_file.split('_')[4]))
        sys.stdout.write(">>> Calculating CLR...\n")
        
        #file initialisation
        all_db = pd.read_csv(args.db, usecols = ['# Qname', 'species', 'genus'], sep='\t')
        mapstat_df = pd.read_csv(mapstat_file, sep='\t', skiprows=(6), usecols = ['# refSequence','fragmentCountAln'])

        #change names of columns for dfs
        mapstat_df = mapstat_df.rename(columns={'# id':'# shortID'})
        mapstat_df = mapstat_df.rename(columns={'# refSequence':'# id'})
        all_db = all_db.rename(columns={'# Qname':'# shortID'})
        
        #splits ids for mapstat and res dataframes
        split_ids(mapstat_df)

        #merge dataframes to create one large
        new_df = pd.merge(mapstat_df, all_db, on="# shortID")
        new_df.dropna(subset=['fragmentCountAln'])
        
        #determine taxonomy level
        if args.tax == 'species':
            species_list = new_df['species'].unique().tolist()
            flag = 'species'
            
            for i in range(0,len(new_df)):
                sum_fragment_count(new_df.iloc[i], species_list, species_dict, flag)
            
            #calculate clr
            read_id = mapstat_file.split('_')[4]
            
            #calculate geometric mean
            species_dict = {k: v for k, v in species_dict.items() if v != 0.0}
            vals_list = species_dict.values()
            species_names = list(species_dict.keys())
            vals_list = [float(i) for i in vals_list]
            geom_mean = geometric_mean(vals_list)
            
            for i in range(0, len(vals_list)):
                clr_val = calc_clr(vals_list[i], geom_mean)
                key = species_names[i]
                clr_dict[key] = clr_val
            

        elif args.tax == 'genus':
            genus_list = new_df['genus'].unique().tolist()
            flag = 'genus'
            
            for i in range(0,len(new_df)):
                sum_fragment_count(new_df.iloc[i], genus_list, genus_dict, flag)
            
            #calculate fpkm
            read_id = mapstat_file.split('_')[4]
            
            species_names = species_dict.keys()
            species_vals = species_dict.values()
            
            clr_vals = clr(species_vals)
            
            zip_iterator = zip(species_names, clr_vals)
            clr_dict = dict(zip_iterator)
        
        with open('{}_clr_values_{}.tab'.format(mapstat_file.split('_')[4],flag),'w') as outfile:
            outfile.write('#name\t#clr_value\n')
            for key, value in clr_dict.items():
                outfile.write('{0}\t{1}\n'.format(key,value))
        
        
        
        
     
                
                
                
                
                
                
                
                
                
                
                
                
                
