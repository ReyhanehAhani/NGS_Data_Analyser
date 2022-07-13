import pandas as pd
import argparse
from typing import List
import io

parser = argparse.ArgumentParser(description='Process csv gene files')
parser.add_argument('mother', type=str, help='path for mother csv file')
parser.add_argument('father', type=str, help='path for father csv file')
parser.add_argument('mother_path', type=str, help='path for mother csv file')
parser.add_argument('father_path', type=str, help='path for father csv file')
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
def read_and_filter_dataset(path: str, parent: str):
    df = pd.read_csv(path)
    return filter_using_parent(df, parent)


def read_and_filter_path_dataset(path: str, parent: str) -> pd.DataFrame:
    return filter_path_using_parent(pd.read_csv(path), parent)

def has_father_mother(row): 
    row = row.sort_values('Parent')    
    parents = tuple(row['Parent'].values)
    return ('mother' in parents) and ('father' in parents)

def main():
    # Read datasets
    df_mother = read_and_filter_dataset(args.mother, 'mother')
    df_father = read_and_filter_dataset(args.father, 'father')

    df_mother_path = read_and_filter_path_dataset(args.mother_path, 'mother')
    df_father_path = read_and_filter_path_dataset(args.father_path, 'father')

    # Merge datasets
    # Drop indexes to avoid incorrect indexes
    df_father = df_father.astype(str)
    df_mother = df_mother.astype(str)
    
    df_mother.reset_index(drop=True, inplace=True)
    df_father.reset_index(drop=True, inplace=True)
    df_mother_path.reset_index(drop=True, inplace=True)
    df_father_path.reset_index(drop=True, inplace=True)

    # Concactinate datasets
    df_merged = pd.concat([df_mother, df_father])
    df_path_merged = pd.concat([df_mother_path, df_father_path])

    # Sort by Gene.refGene
    df = df_merged.sort_values('Gene.refGene')
    df_path = df_path_merged.sort_values('Gene.refGene')

    # Columns that need to be matched
    match_columns = ["Het Iranome", "Hom Iranome", "Het Our DB", "Chr", "Start",
                     "End", "Ref", "Alt", "Zygosity", "Gene.refGene", "ExonicFunc.refGene"]

    gene_exceptions = ['frameshift insertion',
                    'frameshift deletion', 'stopgain', 'stoploss', 'splicing']

    # Couple output
    couple = df.groupby(match_columns)[df.columns].filter(has_father_mother)

    # Compound gene
    compound = df[(df.duplicated(subset=['Gene.refGene'], keep=False))]

    # Dangerous gene
    dangerous = df[(df['ExonicFunc.refGene'].isin(gene_exceptions))
                   | (df['ExonicFunc.ensGene'].isin(gene_exceptions))
                   | (df['ExonicFunc.knownGene'].isin(gene_exceptions))
                   | (df['Func.refGene'].isin(gene_exceptions))
                   | (df['Function_description'].isin(gene_exceptions))]

    # Couple path
    couple_path = df_path.groupby(match_columns)[df_path.columns].filter(has_father_mother)

    # Not shared path
    not_shared_path = df_path[~(df_path.duplicated(
        subset=match_columns, keep=False))]

    # For check in mother
    mother_hom_check = df_mother_path[(df_mother_path['Zygosity'] == 'hom')]

    # For check in father
    father_hom_check = df_father_path[(df_father_path['Zygosity'] == 'hom')]

    print('Writing xlsx ...')

    current_row = 1
    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')

    couple.to_excel(writer, sheet_name='Sheet1', index=False,
                    startrow=current_row, startcol=0)
    current_row += couple.shape[0]

    couple_path.to_excel(writer, sheet_name='Sheet1', header=False,
                         index=False, startrow=current_row, startcol=0)
    current_row += couple_path.shape[0]

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

    worksheet.merge_range(0, 0, 0, 3, 'موارد مشترک در زوج', merge_format)
    current_row += 1

    worksheet.merge_range(current_row, 0, current_row,
                          3, 'ژن مشترک برای احتمال کامپوند', merge_format)
    current_row += 1
    compound.to_excel(writer, sheet_name='Sheet1', header=False,
                      index=False, startrow=current_row)
    current_row += compound.shape[0]
    current_row += 1

    worksheet.merge_range(current_row, 0, current_row,
                          3, 'موارد خطرناک در هر یک از زوجین', merge_format)
    current_row += 1
    dangerous.to_excel(writer, sheet_name='Sheet1', header=False,
                       index=False, startrow=current_row)
    current_row += dangerous.shape[0]
    current_row += 1

    worksheet.merge_range(current_row, 0, current_row,
                          3, 'موارد پاتوژن غیرمشترک', merge_format)
    current_row += 1
    not_shared_path.to_excel(writer, sheet_name='Sheet1',  header=False,
                             index=False, startrow=current_row)
    current_row += not_shared_path.shape[0]
    current_row += 1

    worksheet.merge_range(current_row, 0, current_row,
                          3, 'برای بررسی در پدر', merge_format)
    current_row += 1
    father_hom_check.to_excel(
        writer, sheet_name='Sheet1', index=False,  header=False, startrow=current_row)
    current_row += father_hom_check.shape[0]
    current_row += 1

    worksheet.merge_range(current_row, 0, current_row,
                          3, 'برای بررسی در مادر', merge_format)
    current_row += 1
    mother_hom_check.to_excel(
        writer, sheet_name='Sheet1', index=False,  header=False, startrow=current_row)
    current_row += mother_hom_check.shape[0]
    current_row += 1

    writer.save()


if __name__ == "__main__":
    main()
