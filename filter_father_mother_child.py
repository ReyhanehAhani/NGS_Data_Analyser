import pandas as pd
import argparse
from typing import List
import io

parser = argparse.ArgumentParser(description='Process csv gene files')
parser.add_argument('mother', type=str, help='path for mother csv file',
                    const=1, nargs='?', default='father_mother_child/filtered_MOT18842_2.csv')
parser.add_argument('father', type=str, help='path for mother csv file',
                    const=1, nargs='?', default='father_mother_child/filtered_MOT18839_2.csv')
parser.add_argument('child', type=str, help='path for child csv file',
                    const=1, nargs='?', default='father_mother_child/filtered_MOT18847_3.csv')
parser.add_argument('output', type=str, help='path for output xlsx file',
                    const=1, nargs='?', default='father_mother_child/output.xlsx')
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


# Read csv file and filter using `filter_using_parent`
def read_and_filter_dataset(path: str, parent: str) -> pd.DataFrame:
    return filter_using_parent(pd.read_csv(path), parent)

def has_father_mother_child(row): 
    row = row.sort_values('Parent')    
    return (tuple(row['Parent'].values) == ('child', 'father', 'mother')) and (tuple(row['Zygosity'].values) == ('hom', 'het', 'het'))

def main():
    df_child = read_and_filter_dataset(args.child, 'child')
    df_mother = read_and_filter_dataset(args.mother, 'mother')
    df_father = read_and_filter_dataset(args.father, 'father')

    df_mother.reset_index(drop=True, inplace=True)
    df_father.reset_index(drop=True, inplace=True)
    df_child.reset_index(drop=True, inplace=True)

    df = pd.concat([df_father, df_mother, df_child])
    df = df.sort_values(['Gene.refGene'])

    shared_gene = df.groupby('Gene.refGene')[df.columns].filter(has_father_mother_child)
    not_shared_gene = df[(df.duplicated(subset=['Gene.refGene'])) & (df['Parent'] == 'child')]

    print('Writing xlsx ...')
    current_row = 1
    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')

    shared_gene.to_excel(writer, sheet_name='Sheet1', index=False,
                    startrow=current_row, startcol=0)
    current_row += shared_gene.shape[0]

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

    worksheet.merge_range(current_row+1, 0, current_row+1,
                          3, 'موارد غیرمشترک در فرزند (‌فقط فرزند)', merge_format)
    current_row += 1
    not_shared_gene.to_excel(writer, sheet_name='Sheet1',
                      index=False, startrow=current_row+1)
    current_row += not_shared_gene.shape[0]
    current_row += 1

    writer.save()

if __name__ == "__main__":
    main()