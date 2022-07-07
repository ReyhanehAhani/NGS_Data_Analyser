import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Process csv gene files')
parser.add_argument('mother', type=str, help='path for mother csv file', const=1, nargs='?', default='father_mother/filtered_MOT_T_4904_2.csv')
parser.add_argument('father', type=str, help='path for father csv file', const=1, nargs='?', default='father_mother/filtered_MOT_T_4905_2.csv')
parser.add_argument('output', type=str, help='path for output csv file', const=1, nargs='?', default='father_mother/output.csv')
args = parser.parse_args()


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


# Read datasets
df_mother = filter_using_parent(pd.read_csv(
    args.mother, on_bad_lines='skip'), 'mother')
df_father = filter_using_parent(pd.read_csv(
    args.father, on_bad_lines='skip'), 'father')

# Merge datasets
# Drop indexes to avoid incorrect indexes
df_mother.reset_index(drop=True, inplace=True)
df_father.reset_index(drop=True, inplace=True)
# Concactinate datasets
df_merged = pd.concat([df_mother, df_father])

# Sort by Gene.refGene
df = df_merged.sort_values('Gene.refGene')

df = df.loc[:, df.columns != 'Parent']