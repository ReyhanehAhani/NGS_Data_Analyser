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
    parents = row['Parent'].values
    return (('child' in parents) and ('mother' in parents) and ('father' in parents) and (row['Parent'].size == 3))

def main():
    df_child = read_and_filter_dataset(args.child, 'child')
    df_mother = read_and_filter_dataset(args.mother, 'mother')
    df_father = read_and_filter_dataset(args.father, 'father')

    df_mother.reset_index(drop=True, inplace=True)
    df_father.reset_index(drop=True, inplace=True)
    df_child.reset_index(drop=True, inplace=True)

    df = pd.concat([df_father, df_mother, df_child])
    df = df.sort_values(['Gene.refGene'])

    df.groupby('Gene.refGene')[df.columns].filter(has_father_mother_child).to_csv('output.csv')
    

if __name__ == "__main__":
    main()
