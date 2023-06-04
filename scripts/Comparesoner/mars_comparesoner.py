import pandas as pd
import argparse
import os
import sys
from Bio import SeqIO
import subprocess
from tqdm import tqdm
from datetime import datetime
#bulk variant working with whole folder

parser = argparse.ArgumentParser()

args = parser.add_argument('-path', help='Path to blast_results')
args = parser.add_argument('-evalue', default=1e-3, help='Allow to chose e-value, default = 1e-3')
args = parser.add_argument('-db', default='/home/aster/db/swissprot_16.05.23/uniprot_sprot.fasta')
args = parser.add_argument('-concat', default='/home/aster/Project_Geoarchaeota_28.11/results/mars_comarison/blast_db_concat_Mars_all/mars_concatenated.faa')
args = parser.add_argument('-rm', default=None, help='Allow to choose delite or not bash script for blast, defoult None' )
args = parser.add_argument('-word_size', default=3, help='Allow to choose word size for blast' )

args = parser.parse_args()


def blast_res_reader(path_to_table:str) -> pd.DataFrame:
    try:
        data = pd.read_csv(path_to_table, 
                        sep='\t', header=None)
        data = data.sort_values(by=6)
        data  = data.drop_duplicates(subset=1)
        return data
    except pd.errors.EmptyDataError:
        return False

def swiss_blast_res_reader(path_to_table:str)-> pd.DataFrame:
    try:
        data = pd.read_csv(path_to_table, 
                        sep='\t', header=None)
        data = data.sort_values(by=6)
        data  = data.drop_duplicates(subset=0)
        return data
    
    except pd.errors.EmptyDataError:
        return False


def make_swiss_blastp_script(evalue:float, script_path:str, query: str,
                res_path: str):
    
    with open(script_path, 'w') as write_file:
        write_file.write(f"""\
        #!/bin/bash
        echo 'making blast of -- {query}'
        blastp -query {query} -db /home/aster/db/swissprot_16.05.23/uniprot_sprot.fasta -num_threads 16 -evalue {evalue} -word_size {args.word_size} -outfmt "6 qaccver stitle pident length mismatch gapopen evalue qcovs bitscore" -out {res_path}
        """)
    subprocess.run(['chmod', '+x', script_fpath])  # since script is unexecutable after creating
    subprocess.run([script_fpath], shell=True)

def make_fasta(table:pd.DataFrame, 
               path_to_concat:str, 
               output_path:str, 
               filename:str):

    with open(os.path.join(output_path, f'{filename}.fasta'), 'a') as write_file:
            
        for seq in SeqIO.parse(open(path_to_concat), 'fasta'):
            for hit in table.itertuples():
                if hit[2] in seq.name:
                    write_file.write(f">{hit[2]}\n{seq.seq}\n")


if __name__ == '__main__':

    start = datetime.now()

    swiss_folder = os.path.join('/'.join(str(args.path).split('/')[:-1]), 'swiss_blast') # make separate dir for swissprot blast res
    swiss_res_path = os.path.join(swiss_folder, 'blast_siwss_res')
    swiss_fasta_path = os.path.join(swiss_folder, 'fasta_to_blast')

    fasta_check = False

    ## If folder alrady exixt it fill skip fasta creation since it O(n+nk)
    try:
        os.mkdir(swiss_folder)
    except FileExistsError:
        print('folder already exist')
    try:
        os.mkdir(swiss_res_path)
    except FileExistsError:
        print('folder already exist')
    try:
        os.mkdir(swiss_fasta_path)
    except FileExistsError:
        print('folder already exist')
        fasta_check = True

    for name in tqdm(os.listdir(args.path), desc='Blast filtration and fasta creation'): #iteration throught blasy hits
        if fasta_check:
            break

        fpath = os.path.join(args.path, name)
        table_data = blast_res_reader(fpath)
        if table_data is False:
            print(f'Emty blast result in {fpath}\n')
            continue
        make_fasta(table = table_data, 
                   output_path = swiss_fasta_path, 
                   filename= name, 
                   path_to_concat = args.concat)


    for fastas in tqdm(os.listdir(swiss_fasta_path), desc='blast to swiss prot database'):

        script_fpath = os.path.join(swiss_res_path, f'{fastas}.sh')
        res_path = os.path.join(swiss_res_path, fastas)
        fasta_fpath = os.path.join(swiss_fasta_path, fastas)

        make_swiss_blastp_script(evalue=args.evalue,
            res_path = f'{res_path}.swiss_res.tsv',
            script_path=script_fpath,
            query=fasta_fpath)

        if args.rm:
            os.remove(script_fpath)


    for name in tqdm(os.listdir(swiss_res_path), desc='swiss blast results filtration'):
        if 'sh' not in name:
            fpath = os.path.join(swiss_res_path, name)
            table_data = swiss_blast_res_reader(fpath)
            if table_data is False:
                print(f'Emty blast result in {fpath}\n')
                continue
            
            else:

                table_data.to_csv(fpath, sep='\t', header=None, index=False)

    print(f'Script working time is: {datetime.now() - start}')