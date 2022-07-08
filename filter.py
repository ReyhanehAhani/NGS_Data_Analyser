import pandas as pd
import argparse
from typing import List
import io

parser = argparse.ArgumentParser(description='Process csv gene files')
parser.add_argument('mother', type=str, help='path for mother csv file',
                    const=1, nargs='?', default='father_mother/filtered_MOT_T_4904_2.csv')
parser.add_argument('father', type=str, help='path for father csv file',
                    const=1, nargs='?', default='father_mother/filtered_MOT_T_4905_2.csv')
parser.add_argument('mother_path', type=str, help='path for mother csv file',
                    const=1, nargs='?', default='father_mother/filtered_MOT_T_4904_path_3.csv')
parser.add_argument('father_path', type=str, help='path for father csv file',
                    const=1, nargs='?', default='father_mother/filtered_MOT_T_4905_path_3.csv')
parser.add_argument('output', type=str, help='path for output xlsx file',
                    const=1, nargs='?', default='father_mother/output.xlsx')
args = parser.parse_args()

# Generate column names dynamicly based on first row, some faulty files have excess columns


def generate_column_names(path: str, delimiter: str = ',') -> List[str]:
    with io.open(path, encoding='utf-8') as f:
        columns = f.readline().strip()
        # Retrun all non empty unique columns
        columns = list(dict.fromkeys(
            [x for x in columns.split(delimiter) if x]))
        return columns


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

    return df



# Read csv file and filter using  `filter_using_parent`
def read_and_filter_dataset(path: str, parent: str) -> pd.DataFrame:
    columns = generate_column_names(path)
    # if any excess columns occur, get the columns til the defiend amount
    return filter_using_parent(pd.read_csv(path, names=columns, engine="python", on_bad_lines=lambda x: x[:len(columns)]), parent)

def read_and_filter_path_dataset(path: str, parent: str) -> pd.DataFrame:
    columns = generate_column_names(path)
    # if any excess columns occur, get the columns til the defiend amount
    return filter_path_using_parent(pd.read_csv(path, engine="python"), parent)

def main():
    # Read datasets
    df_mother = read_and_filter_dataset(args.mother, 'mother')
    df_father = read_and_filter_dataset(args.father, 'father')
    df_mother_path = read_and_filter_path_dataset(args.mother_path, 'mother')
    df_father_path = read_and_filter_path_dataset(args.father_path, 'father')
    
    # Merge datasets
    # Drop indexes to avoid incorrect indexes
    df_mother.reset_index(drop=True, inplace=True)
    df_father.reset_index(drop=True, inplace=True)
    df_mother_path.reset_index(drop=True, inplace=True)
    df_father_path.reset_index(drop=True, inplace=True)
    
    # Concactinate datasets
    df_merged = pd.concat([df_mother, df_father])

    # Sort by Gene.refGene
    df = df_merged.sort_values('Gene.refGene')

    # Columns that need to be matched
    match_columns = ["Het Iranome", "Hom Iranome", "Het Our DB", "Chr", "Start",
                     "End", "Ref", "Alt", "Zygosity", "Gene.refGene", "ExonicFunc.refGene"]

    gene_exceptions = ['frameshift insertion',
                       'frameshift deletion', 'stopgain', 'stoploss', 'splice']
    exception_columns = ['ExonicFunc.refGene', 'ExonicFunc.ensGene',
                         'ExonicFunc.knownGene', 'Func.refGene', 'Function_description']

    # Couple output
    couple = df[((df.duplicated(subset=match_columns, keep=False))
                 & (df['Parent'] != df['Parent'].shift(1)))]

    # compound gene
    compound = df[(df.duplicated(subset=['Gene.refGene'], keep=False))]

    # Dangerous gene
    dangerous = df[(df[exception_columns].isin(gene_exceptions))]

    print('Writing xlsx file ...')

    writer = pd.ExcelWriter(args.output, engine='xlsxwriter')
    couple.to_excel(writer, sheet_name='Sheet1',
                    index=False, startcol=0, startrow=1)
    worksheet = writer.sheets['Sheet1']
    workbook  = writer.book

    merge_format = workbook.add_format({
        'bold':     True,
        'align':    'center',
        'valign':   'vcenter',
        'fg_color': '#D7E4BC',
    })

    worksheet.merge_range(0, 0, 0, 3, 'موارد مشترک در زوج', merge_format)

    compound.to_excel(
        writer, sheet_name='Sheet1', index=False, startcol=0, startrow=couple.shape[0]+3)

    worksheet.merge_range(couple.shape[0]+2, 0, couple.shape[0]+2, 3, 'ژن مشترک برای احتمال کامپوند', merge_format)
    
    dangerous.to_excel(
        writer, sheet_name='Sheet1', index=False, startcol=0, startrow=compound.shape[0]+4) #موارد خطرناک در هر یک از زوجین

    worksheet.merge_range(compound.shape[0]+3, 0, compound.shape[0]+3, 3, 'ژن مشترک برای احتمال کامپوند', merge_format)
    

    writer.save()


if __name__ == "__main__":
    main()
