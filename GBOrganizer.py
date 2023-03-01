#!/usr/bin/env python3.10

#Last Updated: 23 Feb 2023
#Author: Pedro Mendes de Souza; pedromsouza0@gmail.com
#Usage: python GBOrganizer.py
#Options: python GBOrganizer.py --help

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from argparse import RawTextHelpFormatter,SUPPRESS
import argparse, re, os, sys
import pandas as pd


# Guide Msg

def guide_msg():
    return "\n\nExemple of usage: pyhton GBOrganizer.py -i sequence.gb -sf Yes -l 1000"

# Author info

def author_info(param):

	valid = 0

	author = ('\n\n\t Any comments or questions? Please email Pedro Souza (author) at'\
	' pedromsouza0@gmail.com\n\n')

	if param.author == True:
		print (author)
		valid += 1
	
	return valid

# Parses and Check 

def get_parameters():

    parser = argparse.ArgumentParser(description=
	'\n\n\nThis script will colect information from GenBank file and organize the data in a excel dataset, and if need divide gene sequences in diferente files.\n'+guide_msg(), 
	usage=SUPPRESS,formatter_class=RawTextHelpFormatter)

    required_param_group = parser.add_argument_group('Required Options')

    required_param_group.add_argument('--input_file','-i', action='store',
	help="\n GenBank file (.gb).\n\n")

    optional_arg_group = parser.add_argument_group('Options')

    optional_arg_group.add_argument('-author', action='store_true',
	help=' \n Print author contact information \n\n')

    optional_arg_group.add_argument('--split_file','-sf', action='store', default= "No", help='\nPut Yes if you want to divide different genes sequences into separate files.\n')

    optional_arg_group.add_argument('--lenght', '-l', action='store', default='0', help='\nUser-difined minimum sequence lenght')


    if len(sys.argv[1:]) == 0:
        print (parser.description)
        print ('\n')
        sys.exit()

    param = parser.parse_args()
	
    quit_eval = author_info(param)
    if quit_eval > 0:
        sys.exit()

    param = parser.parse_args()

    return param

# Main

def GBOrganizer(param):
    data = list(SeqIO.parse(param.input_file, 'gb'))
    data_dict = SeqIO.to_dict(SeqIO.parse(param.input_file, 'gb'))

    # Table

    uni =[]
    more = []

    for i in data:
        if i.features[0].location != i.features[1].location:
            more.append(i)
        else:
            uni.append(i)

    acession = []
    organism = []
    product = []
    leng = []
    journal = []
    authors = []
    location = []

    for i in data:
        for x in i.features:
            if x.qualifiers.get('product'):
                acession.append(i.id)
                product.append(str(x.qualifiers.get('product')).replace("['",'').replace("']",''))
                leng.append(str(x.location).replace("](+)",'').replace('[',"").replace("<",'').replace('>','').replace('](-)','').replace('join{','').replace('}',''))
                organism.append(i.annotations['organism'])
                journal.append(i.annotations['references'][0].journal)
                authors.append(i.annotations['references'][0].authors)
                location.append(str(x.location).replace("](+)",'').replace('[',"").replace("<",'').replace('>','').replace('](-)','').replace('join{','').replace('}',''))

    new_nm = []

    for nm in product:
        if re.match('18S .+', nm) or re.match('small subunit ribosomal',nm):
            new_nm.append('small subunit ribosomal RNA')
        elif re.match('28S .+', nm):
            new_nm.append('large subunit ribosomal RNA')
        else:
            new_nm.append(nm)
            # print(nm)
    product = new_nm


    pre_df = [acession, organism, product, leng, journal, authors, location]
    df = pd.DataFrame(pre_df).transpose()
    df.columns = ['Acession', 'Organism', 'Marker', 'Lenght', 'Journal', 'Authors','Location']

    new_len = []

    for i in df.Lenght:
        if i.count(':') > 1:
            all_frag = i.split(',')
            sum_frag = 0

            for frag in all_frag:
                x,y = frag.split(':')
                z = int(y) - int(x)
                sum_frag += z

            new_len.append(sum_frag)
        else:
            a,b = i.split(':')
            c = int(b) - int(a)
            new_len.append(c)

    df.Lenght = new_len

    df = df.drop_duplicates()
    df = df.reset_index()

    exp = int(param.lenght)
    
    if exp != 0:

        len_in = len(df)

        for i in range(len(df)):
            obs = int(df.loc[i,'Lenght'])
        
            if obs < exp:
                df.drop(index=i,inplace=True)
        
        len_out = len(df)

        print(f'\nOf the {len_in} sequence collected from the .gb file only {len_out} have more then {exp} nucleotides.\n')


    df.to_excel('GB_Organized.xlsx')

    if param.split_file == "Yes":
        
        #Fasta
        path = os.getcwd()

        if not os.path.isdir(f"{path}//Fasta"):
            os.makedirs(f"{path}//Fasta")

        markers = df.Marker.unique()
        # print(markers)

        for marker in markers:
            file = open(str('Fasta/'+marker+'.fasta'),'w')
            df_temp = df[df.Marker == marker]
            df_temp = df_temp.reset_index()
            for i in range(len(df_temp)):

                if df_temp.loc[i,'Location'].count(':') == 1: 

                    start, finish =df_temp.loc[i,'Location'].split(':')
                    sqn = data_dict[str(df_temp.loc[i,'Acession'])].seq[int(start):int(finish)]

                else:
                    locs = df_temp.loc[i,'Location'].split(',')
                    sqn = '$'
                    for l in locs:
                        start, finish = l.split(':')
                        sqn += str(data_dict[str(df_temp.loc[i,'Acession'])].seq[int(start):int(finish)])
                    sqn = sqn.replace('$','')

                file.write('>'+df_temp.loc[i,'Acession']+'_'+(df_temp.loc[i,'Organism'].replace(' ','_'))+'\n'+str(sqn)+'\n')
            file.close()
        #     tamanho += int(len(df_temp))


def run():
    param = get_parameters()
    GBOrganizer(param)

run()
