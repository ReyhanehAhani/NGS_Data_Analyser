import pandas as pd
import argparse
from typing import List
import io

parser = argparse.ArgumentParser(description='Process csv gene files')
parser.add_argument('mother', type=str, help='path for mother csv file',
                    const=1, nargs='?', default='father_mother/filtered_MOT_T_4904_2.csv')
parser.add_argument('father', type=str, help='path for father csv file',
                    const=1, nargs='?', default='father_mother/filtered_MOT_T_4905_2.csv')

# TODO: Use better naming for arguemts and files
parser.add_argument('couple', type=str, help='path for couple csv file',
                    const=1, nargs='?', default='father_mother/couple.csv')
parser.add_argument('compount_gene', type=str, help='path for compount_gene csv file',
                    const=1, nargs='?', default='father_mother/compount_gene.csv')
parser.add_argument('dangerous_gene', type=str, help='path for dangerous_gene csv file',
                    const=1, nargs='?', default='father_mother/dangerous_gene.csv')
args = parser.parse_args()

# Generate column names dynamicly based on first row, some faulty files have excess columns
def generate_column_names(path: str, delimiter: str = ',') -> List[str]:
    with io.open(path, encoding='utf-8') as f:
        columns = f.readline().strip()
        # Retrun all non empty unique columns
        columns = list(dict.fromkeys(
            [x for x in columns.split(delimiter) if x]))
        return columns


def filter_using_parent(df: pd.DataFrame, parent: str):
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

    # Drop columns with all NaN
    df = df.dropna(axis=1, how='all')

    return df


# Read csv file and filter using `filter_using_parent`
def read_and_filter_dataset(path: str, parent: str):
    columns = generate_column_names(args.mother)
    # if any excess columns occur, get the columns til the defiend amount
    return filter_using_parent(pd.read_csv(args.mother, names=columns, engine="python", on_bad_lines=lambda x: x[:len(columns)]), parent)


def main():
    # Read datasets
    df_mother = read_and_filter_dataset(args.mother, 'mother')
    df_father = read_and_filter_dataset(args.father, 'father')

    # Merge datasets
    # Drop indexes to avoid incorrect indexes
    df_mother.reset_index(drop=True, inplace=True)
    df_father.reset_index(drop=True, inplace=True)
    # Concactinate datasets
    df_merged = pd.concat([df_mother, df_father])

    # Sort by Gene.refGene
    df = df_merged.sort_values('Gene.refGene')

    # Columns that need to be matched
    match_columns = ["Het Iranome", "Hom Iranome", "Het Our DB", "Chr", "Start",
                     "End", "Ref", "Alt", "Zygosity", "Gene.refGene", "ExonicFunc.refGene"]
    
    # Exonic values that can't be deleted
    gene_exceptions = ['frameshift insertion', 'frameshift deletion', 'stopgain', 'stoploss', 'splice']
    exception_columns = ['ExonicFunc.refGene', 'ExonicFunc.ensGene', 'ExonicFunc.knownGene', 'Func.refGene']

    # Couple output
    df[((df.duplicated(subset=match_columns, keep=False)) & (df['Parent'] != df['Parent'].shift(1)))].to_csv(args.couple, index=False)

    # Compount gene
    df[(df.duplicated(subset=['Gene.refGene'], keep=False))].to_csv(args.compount_gene, index=False)

    # Dangerous gene
    df[(df[exception_columns].isin(gene_exceptions))].to_csv(args.dangerous_gene, index=False)
    

if __name__ == "__main__":
    main()
