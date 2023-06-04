import pandas as pd
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

#TODO: Make something with colors

parser = argparse.ArgumentParser()

parser.add_argument('-dbsub', help='Path to dbsun tsv file from dbCAN')
parser.add_argument('-overview', help='Path to overview tsv file from dbCAN')
parser.add_argument('-dpi', type=int,default=300, help='Graph dpi')
parser.add_argument('-widh', type=int,default=10, help='Graph width')
parser.add_argument('-high', type=int,default=10, help='Graph hight')
parser.add_argument('-genome_name', default='Genome_1', help='Enter name for the genome; Use in plot title and file path')
parser.add_argument('-verbose', default=None, help='Simplyfy CAZy familis')
args = parser.parse_args()


def tab_reader(table_dbsub, tab_save_path, table_overview = None):

    dbsub_parsed = pd.read_csv(f'{table_dbsub}', sep='\t')
    dbsub_parsed = dbsub_parsed[['dbCAN subfam', 'Subfam Composition']]
    if args.verbose is None:
        for idx, row in dbsub_parsed.iterrows():

            parsed_type = [i for i in list(row['dbCAN subfam']) if i.isalpha() and i != 'e']

            dbsub_parsed.loc[idx, 'CAZ_type'] = ''.join(parsed_type)
        
        # Count Cazyms
        count_CAZ_type_set = set()
        caz_lst = dbsub_parsed['CAZ_type'].to_list()
        for i in caz_lst:
            count_CAZ_type_set.add((i, caz_lst.count(i)))

        #Coun different cazyms
        count_CAZ_all_set = set()
        caz_lst = dbsub_parsed['dbCAN subfam'].to_list()
        for i in caz_lst:
            count_CAZ_all_set.add((i, caz_lst.count(i)))

        count_CAZ_type_set = sorted(list(count_CAZ_type_set), key=lambda x: x[0]) # Сортировка по именам, по идеи этоь делает все это робастным
        count_CAZ_all_set = sorted(list(count_CAZ_all_set), key=lambda x: x[0])

        #Save_tabs
        df_CAZ_type_stat = pd.DataFrame({'CAZ':[i[0] for i in count_CAZ_type_set],
                                        'Count': [i[1] for i in count_CAZ_type_set]})
        df_CAZ_all_stat = pd.DataFrame({'CAZ':[i[0] for i in count_CAZ_all_set],
                                        'Count': [i[1] for i in count_CAZ_all_set]})
        
        df = pd.concat([df_CAZ_type_stat, df_CAZ_all_stat])
        df.to_csv(f'{tab_save_path}', index=False)

        return count_CAZ_type_set, count_CAZ_all_set

    else:
        for idx, row in dbsub_parsed.iterrows():

            parsed_type = [i for i in list(row['dbCAN subfam']) if i.isalpha() and i != 'e']

            dbsub_parsed.loc[idx, 'CAZ_type'] = ''.join(parsed_type)
        
        # Count Cazyms
        count_CAZ_type_set = set()
        caz_lst = dbsub_parsed['CAZ_type'].to_list()
        for i in caz_lst:
            count_CAZ_type_set.add((i, caz_lst.count(i)))

        #Coun different cazyms
        count_CAZ_all_set = set()
        dbsub_parsed['dbCAN subfam'] = dbsub_parsed['dbCAN subfam'].apply(lambda x: x.split('_')[0])
        caz_lst = dbsub_parsed['dbCAN subfam'].to_list()
        for i in caz_lst:
            count_CAZ_all_set.add((i, caz_lst.count(i)))

        count_CAZ_type_set = sorted(list(count_CAZ_type_set), key=lambda x: x[0]) # Сортировка по именам, по идеи этоь делает все это робастным
        count_CAZ_all_set = sorted(list(count_CAZ_all_set), key=lambda x: x[0])

        return count_CAZ_type_set, count_CAZ_all_set

def pie_drawer(list_types:list[tuple], list_all:list[tuple],
               path, genome_name:str=args.genome_name, 
               w:int=args.widh, h:int=args.high, dpi:int=args.dpi):
    
    plt.figure(figsize=(w,h), dpi=dpi)

    plt.pie([i[1] for i in list_types], labels=[f'{i[0]}({i[1]})' for i in list_types], 
            radius=0.5, wedgeprops=dict(width=0.5, edgecolor='w'), labeldistance=0.5, textprops={'fontsize': 10}, rotatelabels=270)

    plt.pie([i[1] for i in list_all], labels=[f'{i[0]} ({i[1]})' for i in list_all], 
            radius=1, wedgeprops=dict(width=0.5, edgecolor='w'),rotatelabels=270, labeldistance=0.6, textprops={'fontsize': 9} )

    plt.title(f'Carboxy active enzymes pie chart for {genome_name}', fontsize=20)
    plt.tight_layout()
    plt.savefig(f'{path}/{genome_name}.svg', format = 'svg', dpi=300);
    
def hist_draw():
    pass   



if __name__ == "__main__":
    
    out_img_folder = Path('./CAZ_sat_output/img') #Create folder for img output data
    out_text_folder = Path('./CAZ_sat_output/text_files') #Create folder for txt output data
    out_img_folder.mkdir(parents=True, exist_ok=True)
    out_text_folder.mkdir(parents=True, exist_ok=True)

    count_CAZ_type_set, count_CAZ_all_set  = tab_reader(args.dbsub, f'./CAZ_sat_output/text_files/{args.genome_name}.csv')
    print(count_CAZ_all_set)
    pie_drawer(count_CAZ_type_set, count_CAZ_all_set, './CAZ_sat_output/img')
