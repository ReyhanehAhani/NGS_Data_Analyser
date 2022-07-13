import pandas as pd
import argparse
from typing import List
import io

parser = argparse.ArgumentParser(description='Process csv gene files')
parser.add_argument('mother', type=str, help='path for mother csv file')
parser.add_argument('father', type=str, help='path for father csv file')
parser.add_argument('child', type=str, help='path for child csv file')
parser.add_argument('mother_path', type=str, help='path for mother path csv file')
parser.add_argument('father_path', type=str, help='path for father path csv file')
parser.add_argument('child_path', type=str, help='path for child path csv file')
parser.add_argument('output', type=str, help='path for output xlsx file')
args = parser.parse_args()


def filter_using_parent(df: pd.DataFrame, parent: str) -> pd.DataFrame:
    # Keep Hom Iranome if it's NaN or 0
    df = df[(df['Hom Iranome'] == '0') | (df['Hom Iranome'] == '.')]

    # Keep refGene if isn't intronic
    df = df[df['Func.refGene'] != 'intronic']

    # Convert '.' to NaN, select columns that are NaN or less than 80, convert NaN back to '.'
    df['Het Iranome'] = pd.to_numeric(df['Het Iranome'], errors='coerce')
    df = df[(df['Het Iranome'].isna()) | (df['Het Iranome'] < 80)]
    df['Het Iranome'] = df['Het Iranome'].fillna('.')

    # Convert '.' to NaN, select columns that are NaN or less than 40, convert NaN back to '.'
    df['Het Our DB'] = pd.to_numeric(df['Het Our DB'], errors='coerce')
    df = df[(df['Het Our DB'].isna()) | (df['Het Our DB'] < 40)]
    df['Het Our DB'] = df['Het Our DB'].fillna('.')

    # Insert parent column
    zygosity_index = df.columns.get_loc('Zygosity')
    df.insert(zygosity_index, 'Parent', parent)

    return df

def filter_path_using_parent(df: pd.DataFrame, parent: str) -> pd.DataFrame:
    # Keep Hom Iranome if it's NaN or 0
    df = df[(df['Hom Iranome'] == '0') | (df['Hom Iranome'] == '.')]
    df = df[~(df['CLNSIG'].str.contains('Benign', case=False))]

    df = df[((df['ValueInfo2'].str.split(':').str[0] == '0/1')
             | (df['ValueInfo2'].str.split(':').str[0] == '1/1'))
            & (df['ValueInfo2'].str.split(':').str[1].str.split(',').str[1].astype(int) > 7)]

    # Insert parent column
    zygosity_index = df.columns.get_loc('Zygosity')
    df.insert(zygosity_index, 'Parent', parent)

    return df

# Read csv file and filter using `filter_using_parent`
def read_and_filter_dataset(path: str, parent: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
    except:
        print('Warning: File contains bad lines, Some data might get corrupted!')
        df = pd.read_csv(path, engine='python', on_bad_lines=lambda x: x[:273])
        df = df.shift(2, axis='columns').fillna('.')

    return filter_using_parent(df, parent)

def read_and_filter_path_dataset(path: str, parent: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
    except:
        print('Warning: File contains bad lines, Some data might get corrupted!')
        df = pd.read_csv(path, engine='python', on_bad_lines=lambda x: x[:273])
        df = df.shift(2, axis='columns').fillna('.')

    return filter_path_using_parent(df, parent)

def has_father_mother_child(row): 
    row = row.sort_values('Parent')    
    zygosity = tuple(row['Zygosity'].values)
    parents = tuple(row['Parent'].values)


    if ('child' in parents and (('father' in parents) or ('mother' in parents))):
        if zygosity[0] == 'hom' and zygosity[1] == 'het':
            if len(zygosity) > 2: 
                if zygosity[2] == 'het': 
                    return True
                else:
                    False
            return True
    
    return False


def not_shared_father_mother_child(row): 
    row = row.sort_values('Parent')
    return ('father' not in row['Parent'].values) and ('mother' not in row['Parent'].values)


def main():
    df_child = read_and_filter_dataset(args.child, 'child')
    df_mother = read_and_filter_dataset(args.mother, 'mother')
    df_father = read_and_filter_dataset(args.father, 'father')

    df_child_path = read_and_filter_path_dataset(args.child_path, 'child')
    df_mother_path = read_and_filter_path_dataset(args.mother_path, 'mother')
    df_father_path = read_and_filter_path_dataset(args.father_path, 'father')

    df_mother.reset_index(drop=True, inplace=True)
    df_father.reset_index(drop=True, inplace=True)
    df_child.reset_index(drop=True, inplace=True)

    df_mother_path.reset_index(drop=True, inplace=True)
    df_father_path.reset_index(drop=True, inplace=True)
    df_child_path.reset_index(drop=True, inplace=True)

    df = pd.concat([df_father, df_mother, df_child])
    df = df.sort_values(['Gene.refGene'])

    df_path = pd.concat([df_father_path, df_mother_path, df_child_path])
    df_path = df_path.sort_values(['Gene.refGene'])

    shared_gene = df.groupby('Gene.refGene')[df.columns].filter(has_father_mother_child)
    shared_gene_path = df_path.groupby('Gene.refGene')[df_path.columns].filter(has_father_mother_child)

    not_shared_gene = df.groupby('Gene.refGene')[df.columns].filter(not_shared_father_mother_child)
    not_shared_gene_path = df_path.groupby('Gene.refGene')[df_path.columns].filter(not_shared_father_mother_child)

    print('Writing xlsx ...')
    current_row = 1
    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')

    shared_gene.to_excel(writer, sheet_name='Sheet1', index=False,
                    startrow=current_row, startcol=0)
    current_row += shared_gene.shape[0]

    shared_gene_path.to_excel(writer, sheet_name='Sheet1', index=False, header=False,
                    startrow=current_row+1, startcol=0)
    current_row += shared_gene_path.shape[0]

    worksheet = writer.sheets['Sheet1']
    workbook = writer.book

    merge_format = workbook.add_format({
        'bold':     True,
        'align':    'center',
        'valign':   'vcenter',
        'fg_color': 'black',
        'font_color': 'white',
        'font_size': 16
    })

    worksheet.merge_range(0, 0, 0, 3, 'موارد مشترک در پدر و مادر و فرزند', merge_format)
    current_row += 1

    worksheet.merge_range(current_row, 0, current_row,
                          3, 'موارد غیرمشترک در فرزند (‌فقط فرزند)', merge_format)
    current_row += 1

    not_shared_gene.to_excel(writer, sheet_name='Sheet1',
                      index=False, startrow=current_row,  header=False)
    current_row += not_shared_gene.shape[0]

    not_shared_gene_path.to_excel(writer, sheet_name='Sheet1',
                      index=False, startrow=current_row,  header=False)
    current_row += not_shared_gene_path.shape[0]

    writer.save()

if __name__ == "__main__":
    main()
